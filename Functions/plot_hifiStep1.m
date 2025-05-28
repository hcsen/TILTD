function [c, ceq] = plot_hifiStep1(x,consts)

% Constants
Np = consts(1);         % Number of phases
N = consts(2);          % Number of segments per phase
N_perts = consts(3);            % Number of perturbing bodies
N_ephem = consts(4);            % Number of point in ephemeris data
N_thrust = consts(5);           % Number of phases with thrust in
N_flybys = consts(6);
NpCurrent = consts(7);          % What number phase currently in
mu_central = consts(8);         % Central body standard gravitational parameter in DU^3/TU^2
r_central = consts(9);          % Central body radius in DU
SOI = consts(10);                % Central body sphere of influence radius in DU
m0 = consts(11);                % Initial spacecraft mass in MU
N_act = consts(12);             % Number of active thrusters
D = consts(13);                 % Thrust duty cycle
Tmax = consts(14);              % Max thrust in MU DU/TU^2 (from kN)
Isp = consts(15);               % Specific impulse in TU
g0 = consts(16);                % Acceleration due to gravity on Earth's surface in DU/TU^2
DU = consts(17);                % Distance unit in km
TU = consts(18);                % Time unit in s
MU = consts(19);                % Mass unit in kg
whichThrust = consts(20 : 19+N_thrust);     % Index of phases with thrust
flybySequence = consts(20+N_thrust : 19+N_thrust+N_flybys);   % Sequence of flybys
startBody = consts(20+N_thrust+N_flybys);
endBody = consts(21+N_thrust+N_flybys);
whichFlyby = consts(22+N_thrust+N_flybys : 21+N_thrust+2*N_flybys);
mu_perts = consts(22+N_thrust+2*N_flybys : 21+N_thrust+2*N_flybys+N_perts); % Standard gravitational parameter of perturbing bodies in DU^3/TU^2
r_flybys = consts(22+N_thrust+2*N_flybys+N_perts : 21+N_thrust+3*N_flybys+N_perts);   % Radius of flyby bodies in DU
h_mins = consts(22+N_thrust+3*N_flybys+N_perts : 21+N_thrust+4*N_flybys+N_perts);   % Min flyby heights in DU
h_maxs = consts(22+N_thrust+4*N_flybys+N_perts : 21+N_thrust+5*N_flybys+N_perts);   % Max flyby heights in DU
et = consts(22+N_thrust+5*N_flybys+N_perts : 21+N_thrust+5*N_flybys+N_perts+N_ephem);   % Epochs at which ephemerides taken in DU
stateLong = zeros(6*N_perts, N_ephem);
for i = 1:N_perts
    for j = 1:6
        stateLong(j+(i-1)*6,:) = consts(22+N_thrust+5*N_flybys+N_perts+(j+6*(i-1))*N_ephem : 21+N_thrust+5*N_flybys+N_perts+(j+6*(i-1)+1)*N_ephem);
    end
end

flybyBodies = unique(flybySequence);
N_fbodies = length(flybyBodies);
N_int = consts(22+N_thrust+5*N_flybys+N_perts+(7+6*(N_perts-1))*N_ephem);
SOI_flybys = consts(23+N_thrust+5*N_flybys+N_perts+(7+6*(N_perts-1))*N_ephem : 22+N_thrust+5*N_flybys+N_perts+(7+6*(N_perts-1))*N_ephem+N_fbodies);
DU_flybys = consts(23+N_thrust+5*N_flybys+N_perts+(7+6*(N_perts-1))*N_ephem+N_fbodies : 22+N_thrust+5*N_flybys+N_perts+(7+6*(N_perts-1))*N_ephem+2*N_fbodies);
TU_flybys = consts(23+N_thrust+5*N_flybys+N_perts+(7+6*(N_perts-1))*N_ephem+2*N_fbodies : 22+N_thrust+5*N_flybys+N_perts+(7+6*(N_perts-1))*N_ephem+3*N_fbodies);
MU_flybys = consts(23+N_thrust+5*N_flybys+N_perts+(7+6*(N_perts-1))*N_ephem+3*N_fbodies : 22+N_thrust+5*N_flybys+N_perts+(7+6*(N_perts-1))*N_ephem+4*N_fbodies);

VU = DU/TU;
VU_flybys = DU_flybys./TU_flybys;
thrustPresent = length(find(whichThrust<=NpCurrent));
flybyPresent = find(whichFlyby<=NpCurrent,1,'last');

FBEnd = consts(23+N_thrust+5*N_flybys+N_perts+(7+6*(N_perts-1))*N_ephem+4*N_fbodies : 22+N_thrust+5*N_flybys+N_perts+(7+6*(N_perts-1))*N_ephem+4*N_fbodies+flybyPresent);

%% Decision variables

% For initial phase: initial epoch, initial RA DEC, initial vinf relative to Jupiter
% For each phase: tof, mf, vinfi and vinff for the gravity assist at the
% end of the trajectory, control vector u
% For final phase: vinfi but no vinff because not doing a GA
% Gravity assist phases: initial and final RA DEC, vinfs, rp, vp, tofs

dt = zeros(1, Np);
vinfi = zeros(3, N_flybys);
vinff = zeros(3, N_flybys);
mf = zeros(1, N_thrust);
u = zeros(3*N_thrust, N);

RAi = zeros(1,N_flybys);
DECi = zeros(1,N_flybys);
RAf = zeros(1,N_flybys);
DECf = zeros(1,N_flybys);
rp = zeros(3, N_flybys);
vp = zeros(3, N_flybys);
dt_GAa = zeros(1, N_flybys);
dt_GAb = zeros(1, N_flybys);

%% First phase

t0 = x(1);                                  % Initial epoch in TU
% If trajectory starts on SOI
if startBody == 1
    RA = x(2);                                  % Initial RA in rad
    DEC = x(3);                                 % Initial DEC in rad
    vi = x(4:6);
    ri = radec2cart(SOI, RA, DEC);
    decisionInd = 7;
    
    % Or, if trajectory starts in orbit around central body
elseif startBody == 2
    a = x(2);
    AoP = x(3);
    e = x(4);
    inc = x(5);
    RAAN = x(6);
    t_anom = x(7);
    orbels = [a, AoP, e, inc, RAAN, t_anom];
    [ri, vi] = orb2state(mu_central, orbels);
    decisionInd = 8;
    
    % Or, if trajectory starts at location of perturbing body
else
    ri = spline(et, stateLong((1:3)+((startBody-2)-1)*6, :), t0*TU);
    vi = spline(et, stateLong((4:6)+((startBody-2)-1)*6, :), t0*TU).' + x(2:4);
    decisionInd = 5;
end
thrustInd = 1;
freeInd = 1;
flybyInd = 1;

%% All phases

% If there are non-flyby intermediate control nodes, made some position and
% velocity arrays for them
if N_flybys ~= Np-1
    r_free = zeros(3,Np-1-N_flybys);
    v_free = zeros(3,Np-1-N_flybys);
end

for i = 1:NpCurrent
    dt(i) = x(decisionInd);
    % If phase ends in a flyby
    if any(i == whichFlyby)
        RAi(flybyInd) = x(decisionInd+1);
        DECi(flybyInd) = x(decisionInd+2);
        vinfi(:,flybyInd) = x(decisionInd+3 : decisionInd+5);
        decisionInd = decisionInd+6;
        
        % If it's the final phase
    elseif i == Np
        % If end at SOI
        if endBody == 1
            RAend = x(decisionInd+1);                                  % Initial RA in rad
            DECend = x(decisionInd+2);                                 % Initial DEC in rad
            vf = x(decisionInd+3:decisionInd+5);
            rf = radec2cart(SOI, RAend, DECend);
            decisionInd = decisionInd+6;    
            
            % Or, if trajectory ends in orbit around central body
        elseif endBody == 2
            aend = x(decisionInd+1);
            AoPend = x(decisionInd+2);
            eend = x(decisionInd+3);
            incend = x(decisionInd+4);
            RAANend = x(decisionInd+5);
            t_anomend = x(decisionInd+6);
            orbels = [aend, AoPend, eend, incend, RAANend, t_anomend];
            [rf, vf] = orb2state(mu_central, orbels);
            decisionInd = decisionInd+7;

            % Or, if trajectory ends at location of perturbing body
        else
            rf = spline(et, stateLong((1:3)+((endBody-2)-1)*6, :), (t0+sum(dt))*TU);
            vfinf = x(decisionInd+1:decisionInd+3);
            vf = spline(et, stateLong((4:6)+((endBody-2)-1)*6, :), (t0+sum(dt))*TU).' + vfinf;
            decisionInd = decisionInd+4;
        end
        
        % If it's something else with position and velocity defined
    else
        r_free(:, freeInd) = x(decisionInd+1:decisionInd+3);
        v_free(:,freeInd) = x(decisionInd+4:decisionInd+6);
        decisionInd = decisionInd+7;        
    end
    
    % If this is a thrust arc
    if any(i == whichThrust)
        mf(thrustInd) = x(decisionInd);
        u((1:3)+(thrustInd-1)*3,:) = [x(decisionInd+1 : decisionInd+N); ...
                x(decisionInd+N+1 : decisionInd+2*N); x(decisionInd+2*N+1 : decisionInd + 3*N)];
        decisionInd = decisionInd + 3*N + 1;
        thrustInd = thrustInd + 1;
    end
    
    % If GA is being propagated
    if any(i == whichFlyby) && FBEnd(i) == 1
         dt_GAa(flybyInd) = x(decisionInd);
         dt_GAb(flybyInd) = x(decisionInd+1);
         rp(:,flybyInd) = x(decisionInd+2 : decisionInd+4);
         vp(:,flybyInd) = x(decisionInd+5 : decisionInd+7);
         RAf(flybyInd) = x(decisionInd+8);
         DECf(flybyInd) = x(decisionInd+9);
         vinff(:,flybyInd) = x(decisionInd+10 : decisionInd+12);
         decisionInd = decisionInd+13;
    elseif any(i == whichFlyby)
        vinff(:,flybyInd) = x(decisionInd : decisionInd+2);
        decisionInd = decisionInd+3;
    end

    flybyInd = flybyInd + 1;
end
freePresent = freeInd - 1;

%% Spline for ephemeris data

% Need a position for each perturbing body at the centre of each segment.
% Need velocity of each flyby body at the control node where the GA is
impulse_times = zeros(NpCurrent, N);
pert_allx = zeros(N_perts, NpCurrent*N);
pert_ally = zeros(N_perts, NpCurrent*N);
pert_allz = zeros(N_perts, NpCurrent*N);
flyby_v = zeros(3, N_flybys);
flyby_vStart = zeros(3, N_flybys);
flyby_vEnd = zeros(3, N_flybys);
pert_FBxa = zeros(N_perts, N_flybys*N);
pert_FBya = zeros(N_perts, N_flybys*N);
pert_FBza = zeros(N_perts, N_flybys*N);
pert_FBxb = zeros(N_perts, N_flybys*N);
pert_FByb = zeros(N_perts, N_flybys*N);
pert_FBzb = zeros(N_perts, N_flybys*N);
FB_impulse_timesa = zeros(N_flybys, N*2);
FB_impulse_timesb = zeros(N_flybys, N*2);
stateFBcenta = zeros(3*N_flybys,N*2);
stateFBtoPerta = zeros(3*N_perts*N_flybys,N*2);
stateFBcentb = zeros(3*N_flybys,N*2);
stateFBtoPertb = zeros(3*N_perts*N_flybys,N*2);
mu_FBCent = zeros(1,N_flybys);  % To store central body mu of each flyby
mu_pertFB = zeros(N_perts, N_flybys);       % To store perturbing mus each flyby
flybyInd = 1;
tstart = t0*TU;

for i = 1:NpCurrent
    thalf = (dt(i)/(2*N))*TU;
    if any(whichFlyby == i-1)
        bShort = find(flybyBodies == flybySequence(flybyInd-1));
        tstart = tstart + (dt_GAa(flybyInd-1) + dt_GAb(flybyInd-1))*TU_flybys(bShort);
    end
    tend = tstart + dt(i)*TU;
    impulse_times(i,:) = linspace(tstart+thalf, tend-thalf, N);
    tstart = tend;      % Start of the next phase is end of this phase
    
    % Position values for each perturbing body at each impulse
    for j = 1:N_perts
        pert_allx(j, i*N-(N-1):i*N) = spline(et, stateLong(1+(j-1)*6, :), impulse_times(i,:));
        pert_ally(j, i*N-(N-1):i*N) = spline(et, stateLong(2+(j-1)*6, :), impulse_times(i,:));
        pert_allz(j, i*N-(N-1):i*N) = spline(et, stateLong(3+(j-1)*6, :), impulse_times(i,:));
    end
    
    % If propagating flybys, need to get perturbing bodies from the flyby
    % body's frame
    if any(i == whichFlyby) && FBEnd(i) ==1
        bCurrent = flybySequence(flybyInd);
        bShort = find(flybyBodies == bCurrent);
        % Adding a bit to the times for ephemeris data because sometimes it
        % doesn't cover the whole phase
        thalf = (dt_GAa(flybyInd)/(4*N))*TU_flybys(bShort);
%         FB_impulse_timesa(flybyInd,:) = linspace(tend-thalf, tend + dt_GAa(flybyInd)*TU_flybys(bShort)+thalf, N*2);
        FB_impulse_timesa(flybyInd,:) = linspace(tend-dt_GAa(flybyInd)*TU_flybys(bShort), tend + 2*dt_GAa(flybyInd)*TU_flybys(bShort), N*2);
        thalf = (dt_GAb(flybyInd)/(4*N))*TU_flybys(bShort);
%         FB_impulse_timesb(flybyInd,:) = linspace(tend+dt_GAa(flybyInd)*TU_flybys(bShort)-thalf, ...
%             tend+dt_GAa(flybyInd)*TU_flybys(bShort)+dt_GAb(flybyInd)*TU_flybys(bShort)+thalf, N*2);
        FB_impulse_timesb(flybyInd,:) = linspace(tend+dt_GAa(flybyInd)*TU_flybys(bShort)-dt_GAb(flybyInd)*TU_flybys(bShort), ...
            tend+dt_GAa(flybyInd)*TU_flybys(bShort)+2*dt_GAb(flybyInd)*TU_flybys(bShort), N*2);
        
        % Flyby data relative to main central body
        for j = 1:N_perts
            pert_FBxa(j, flybyInd*N*2-(N*2-1):flybyInd*N*2) = spline(et, stateLong(1+(j-1)*6, :), FB_impulse_timesa(flybyInd,:))*DU/DU_flybys(bShort);
            pert_FBya(j, flybyInd*N*2-(N*2-1):flybyInd*N*2) = spline(et, stateLong(2+(j-1)*6, :), FB_impulse_timesa(flybyInd,:))*DU/DU_flybys(bShort);
            pert_FBza(j, flybyInd*N*2-(N*2-1):flybyInd*N*2) = spline(et, stateLong(3+(j-1)*6, :), FB_impulse_timesa(flybyInd,:))*DU/DU_flybys(bShort);
            pert_FBxb(j, flybyInd*N*2-(N*2-1):flybyInd*N*2) = spline(et, stateLong(1+(j-1)*6, :), FB_impulse_timesb(flybyInd,:))*DU/DU_flybys(bShort);
            pert_FByb(j, flybyInd*N*2-(N*2-1):flybyInd*N*2) = spline(et, stateLong(2+(j-1)*6, :), FB_impulse_timesb(flybyInd,:))*DU/DU_flybys(bShort);
            pert_FBzb(j, flybyInd*N*2-(N*2-1):flybyInd*N*2) = spline(et, stateLong(3+(j-1)*6, :), FB_impulse_timesb(flybyInd,:))*DU/DU_flybys(bShort);
        end
        
        % Main central body as a perturbation relative to flyby
        stateFBcenta(1+(flybyInd-1)*3,:) = -pert_FBxa(bCurrent, flybyInd*N*2-(N*2-1):flybyInd*N*2);
        stateFBcenta(2+(flybyInd-1)*3,:) = -pert_FBya(bCurrent, flybyInd*N*2-(N*2-1):flybyInd*N*2);
        stateFBcenta(3+(flybyInd-1)*3,:) = -pert_FBza(bCurrent, flybyInd*N*2-(N*2-1):flybyInd*N*2);
        stateFBcentb(1+(flybyInd-1)*3,:) = -pert_FBxb(bCurrent, flybyInd*N*2-(N*2-1):flybyInd*N*2);
        stateFBcentb(2+(flybyInd-1)*3,:) = -pert_FByb(bCurrent, flybyInd*N*2-(N*2-1):flybyInd*N*2);
        stateFBcentb(3+(flybyInd-1)*3,:) = -pert_FBzb(bCurrent, flybyInd*N*2-(N*2-1):flybyInd*N*2);
        
        % Now get perturbing bodies from flyby body to them
        for j = 1:N_perts
            if j ~= bCurrent
                stateFBtoPerta(1+(j-1)*3+(flybyInd-1)*N_perts*3,:) = stateFBcenta(1+(flybyInd-1)*3,:) + pert_FBxa(j, flybyInd*N*2-(N*2-1):flybyInd*N*2);
                stateFBtoPerta(2+(j-1)*3+(flybyInd-1)*N_perts*3,:) = stateFBcenta(2+(flybyInd-1)*3,:) + pert_FBya(j, flybyInd*N*2-(N*2-1):flybyInd*N*2);
                stateFBtoPerta(3+(j-1)*3+(flybyInd-1)*N_perts*3,:) = stateFBcenta(3+(flybyInd-1)*3,:) + pert_FBza(j, flybyInd*N*2-(N*2-1):flybyInd*N*2);
                stateFBtoPertb(1+(j-1)*3+(flybyInd-1)*N_perts*3,:) = stateFBcentb(1+(flybyInd-1)*3,:) + pert_FBxb(j, flybyInd*N*2-(N*2-1):flybyInd*N*2);
                stateFBtoPertb(2+(j-1)*3+(flybyInd-1)*N_perts*3,:) = stateFBcentb(2+(flybyInd-1)*3,:) + pert_FByb(j, flybyInd*N*2-(N*2-1):flybyInd*N*2);
                stateFBtoPertb(3+(j-1)*3+(flybyInd-1)*N_perts*3,:) = stateFBcentb(3+(flybyInd-1)*3,:) + pert_FBzb(j, flybyInd*N*2-(N*2-1):flybyInd*N*2);
                mu_pertFB(j,flybyInd) = mu_perts(j)*(DU^3/TU^2)/(DU_flybys(bShort)^3/TU_flybys(bShort)^2);
            end
        end
        
        mu_FBCent(flybyInd) = mu_perts(bCurrent)*(DU^3/TU^2)/(DU_flybys(bShort)^3/TU_flybys(bShort)^2);
        
        % Remove the empty gap where the new central body was in perturbation
        % state
        stateFBtoPertaTemp = nonzeros(stateFBtoPerta((1:3+(N_perts-1)*3)+(flybyInd-1)*N_perts*3,:));
        stateFBtoPertaTemp = reshape(stateFBtoPertaTemp, 3*(N_perts-1),N*2); 
        stateFBtoPerta((1:3+(N_perts-2)*3)+(flybyInd-1)*N_perts*3,:) = stateFBtoPertaTemp;
        stateFBtoPerta((1:3)+(N_perts-1)*3+(flybyInd-1)*N_perts*3,:) = stateFBcenta((1:3)+(flybyInd-1)*3,:);
        stateFBtoPertbTemp = nonzeros(stateFBtoPertb((1:3+(N_perts-1)*3)+(flybyInd-1)*N_perts*3,:));
        stateFBtoPertbTemp = reshape(stateFBtoPertbTemp, 3*(N_perts-1),N*2); 
        stateFBtoPertb((1:3+(N_perts-2)*3)+(flybyInd-1)*N_perts*3,:) = stateFBtoPertbTemp;
        stateFBtoPertb((1:3)+(N_perts-1)*3+(flybyInd-1)*N_perts*3,:) = stateFBcentb((1:3)+(flybyInd-1)*3,:);
        mu_pertFBtemp = nonzeros(mu_pertFB(:,flybyInd)).';
        mu_pertFB(:,flybyInd) = [mu_pertFBtemp, mu_central*(DU^3/TU^2)/(DU_flybys(bShort)^3/TU_flybys(bShort)^2)];
        
        % Velocity of flyby body at flyby epoch, start, and end of velocity
        flyby_v(:,i) = spline(et, stateLong((4:6)+(bCurrent-1)*6,:), tend+dt_GAa(flybyInd)*TU_flybys(bShort));
        flyby_vStart(:,i) = spline(et, stateLong((4:6)+(bCurrent-1)*6,:), tend);
        flyby_vEnd(:,i) = spline(et, stateLong((4:6)+(bCurrent-1)*6,:), tend+(dt_GAa(flybyInd)+dt_GAb(flybyInd))*TU_flybys(bShort));
        
        flybyInd = flybyInd+1;
        
        % If not propagating, still need the flyby velocity
    elseif any(i == whichFlyby)
        flyby_v(:, i) = spline(et, stateLong((4:6)+(flybySequence(i)-1)*6,:), tend);
    end
end

%% Running Sims-Flanagan one phase at a time

x_all = zeros(NpCurrent, N+4);
y_all = zeros(NpCurrent, N+4);
z_all = zeros(NpCurrent, N+4);
vx_all = zeros(NpCurrent, N+4);
vy_all = zeros(NpCurrent, N+4);
vz_all = zeros(NpCurrent, N+4);
m_all = zeros(NpCurrent, N+2);
t_all = zeros(NpCurrent, N+4);
% Don't know how long GA will end up through integrator, keep track
length_GAaF = zeros(1,N_flybys);
length_GAaB = zeros(1,N_flybys);
length_GAbF = zeros(1,N_flybys);
length_GAbB = zeros(1,N_flybys);
delta_all = zeros(1, N_flybys);
rp_all = zeros(1, N_flybys);
sfconstsCoast = [N, mu_central, N_perts, mu_perts];
sfconstsThrust = [N, mu_central, Tmax, Isp, g0, N_act, D, N_perts, mu_perts];
rStart = ri;
vStart = vi;
tStart = t0;
mStart = m0;
mEnd = m0;
thrustInd = 1;
flybyInd = 1;
freeInd = 1;

for i = 1:NpCurrent    
    
    tEnd = tStart + dt(i);
    
    % End values if this is the last phase
    if i == Np
        rEnd = rf;
        vEnd = vf;    
        
    % End values if this is a flyby
    elseif any(i == whichFlyby)
        bCurrent = flybySequence(flybyInd);
        bShort = find(flybyBodies == bCurrent);
        rStartGAa = radec2cart(SOI_flybys(bShort), RAi(i), DECi(i));
        pertxGA = spline(et, stateLong(1+(bCurrent-1)*6, :), tEnd*TU);
        pertyGA = spline(et, stateLong(2+(bCurrent-1)*6, :), tEnd*TU);
        pertzGA = spline(et, stateLong(3+(bCurrent-1)*6, :), tEnd*TU);
        rEnd = rStartGAa*DU_flybys(bShort)/DU + [pertxGA, pertyGA, pertzGA];
        vEnd = flyby_vStart(:,flybyInd) + vinfi(:,flybyInd);
        
    % End values if this is a free point
    else
        rEnd = r_free(:,freeInd);
        vEnd = v_free(:,freeInd);
        freeInd = freeInd + 1;
    end         % End of finding end point
    
%% If this is a thrust arc

    if any(i == whichThrust)
        mEnd = mf(thrustInd);
        
        % Forwards shooting
        [tFwd,rFwd,vFwd,mFwd] = sf_forwards_thrust(rStart, vStart, dt(i), tStart, mStart, u((1:3)+(thrustInd-1)*3,:), ...
            pert_allx(:, i*N-(N-1):i*N), pert_ally(:, i*N-(N-1):i*N), pert_allz(:, i*N-(N-1):i*N), sfconstsThrust);
        x_all(i, 1:N/2+2) = rFwd(:,1);
        y_all(i, 1:N/2+2) = rFwd(:,2);
        z_all(i, 1:N/2+2) = rFwd(:,3);
        vx_all(i, 1:N/2+2) = vFwd(:,1);
        vy_all(i, 1:N/2+2) = vFwd(:,2);
        vz_all(i, 1:N/2+2) = vFwd(:,3);
        m_all(i, 1:N/2+1) = mFwd;
        t_all(i, 1:N/2+2) = tFwd;

        % Backwards shooting
        [tBkd, rBkd, vBkd, mBkd] = sf_backwards_thrust(rEnd, vEnd, dt(i), tEnd, mEnd, u((1:3)+(thrustInd-1)*3,:), ...
            pert_allx(:, i*N-(N-1):i*N), pert_ally(:, i*N-(N-1):i*N), pert_allz(:, i*N-(N-1):i*N), sfconstsThrust);
        x_all(i, N/2+3:end) = rBkd(:,1);
        y_all(i, N/2+3:end) = rBkd(:,2);
        z_all(i, N/2+3:end) = rBkd(:,3);
        vx_all(i, N/2+3:end) = vBkd(:,1);
        vy_all(i, N/2+3:end) = vBkd(:,2);
        vz_all(i, N/2+3:end) = vBkd(:,3);
        m_all(i, N/2+2:end) = mBkd;
        t_all(i, N/2+3:end) = tBkd;
        
        mStart = mEnd;
        thrustInd = thrustInd + 1;

%% If this is a coast arc
    else
        % Forwards shooting
        [tFwd,rFwd,vFwd] = sf_forwards(rStart, vStart, dt(i), tStart, pert_allx(:, i*N-(N-1):i*N), ...
            pert_ally(:, i*N-(N-1):i*N), pert_allz(:, i*N-(N-1):i*N), sfconstsCoast);
        x_all(i, 1:N/2+2) = rFwd(:,1);
        y_all(i, 1:N/2+2) = rFwd(:,2);
        z_all(i, 1:N/2+2) = rFwd(:,3);
        vx_all(i, 1:N/2+2) = vFwd(:,1);
        vy_all(i, 1:N/2+2) = vFwd(:,2);
        vz_all(i, 1:N/2+2) = vFwd(:,3);
        t_all(i, 1:N/2+2) = tFwd;

        % Backwards shooting
        [tBkd,rBkd,vBkd] = sf_backwards(rEnd, vEnd, dt(i), tEnd, pert_allx(:, i*N-(N-1):i*N), ...
            pert_ally(:, i*N-(N-1):i*N), pert_allz(:, i*N-(N-1):i*N), sfconstsCoast);
        x_all(i, N/2+3:end) = rBkd(:,1);
        y_all(i, N/2+3:end) = rBkd(:,2);
        z_all(i, N/2+3:end) = rBkd(:,3);
        vx_all(i, N/2+3:end) = vBkd(:,1);
        vy_all(i, N/2+3:end) = vBkd(:,2);
        vz_all(i, N/2+3:end) = vBkd(:,3);
        t_all(i, N/2+3:end) = tBkd;
    end
    
%% Apply gravity assist, if necessary

    if any(i == whichFlyby)
        vEnd = flyby_v(:,flybyInd) + vinff(:, flybyInd);

        % Turn angle to check it's feasible
        delta_all(flybyInd) = acos((dot(vinfi(:,flybyInd),vinff(:,flybyInd)))/(norm(vinfi(:,flybyInd))*norm(vinff(:,flybyInd))));

        % Flyby periapse radius for constraint on height and turn angle
        rp_all(flybyInd) = mu_perts(flybySequence(flybyInd))/(norm(vinff(:,flybyInd))^2)*(1/(sin(delta_all(flybyInd)/2)) - 1);
        
        % If the gravity assist is being propagated fully, do so with
        % explicit integration
        if FBEnd(flybyInd) == 1
            %% First phase goes from SOI to periapse
            
            vStartGAa = vinfi(:,flybyInd)*VU/VU_flybys(bShort);
            
            x_centPertBig = zeros(N_perts, N*2);
            y_centPertBig = zeros(N_perts, N*2);
            z_centPertBig = zeros(N_perts, N*2);
            for j = 1:N_perts
                x_centPertBig(j, :) = stateFBtoPerta(1+(j-1)*3+(flybyInd-1)*N_perts*3,:);
                y_centPertBig(j, :) = stateFBtoPerta(2+(j-1)*3+(flybyInd-1)*N_perts*3,:);
                z_centPertBig(j, :) = stateFBtoPerta(3+(j-1)*3+(flybyInd-1)*N_perts*3,:);
            end

            % Forwards shooting
            X0 = [rStartGAa, vStartGAa.'];
            
            tspan = [tEnd*TU/TU_flybys(bShort), tEnd*TU/TU_flybys(bShort) + dt_GAa(flybyInd)/2];
            
            [ta,state] = ode45(@(t,X0) EoM_hifi_coast(t, X0, mu_FBCent(flybyInd), N_perts, mu_pertFB(:,flybyInd), x_centPertBig, ...
                y_centPertBig, z_centPertBig, FB_impulse_timesa(i,:), TU_flybys(bShort)), tspan, X0);
            
            % Figure out how many data points are in phase now
            length_GAaF(flybyInd) = length(state(:,1));
            
            % Sort out any difference in size between different vectors
            if flybyInd == 1
                GAa_xF = zeros(N_flybys,length_GAaF(1));
                GAa_yF = zeros(N_flybys,length_GAaF(1));
                GAa_zF = zeros(N_flybys,length_GAaF(1));
                GAa_vxF = zeros(N_flybys,length_GAaF(1));
                GAa_vyF = zeros(N_flybys,length_GAaF(1));
                GAa_vzF = zeros(N_flybys,length_GAaF(1));
                GAa_tF = zeros(N_flybys,length_GAaF(1));
                
            elseif length_GAaF(flybyInd) > maxGAaF
                xTemp = zeros(N_flybys, length_GAaF(flybyInd));
                yTemp = zeros(N_flybys, length_GAaF(flybyInd));
                zTemp = zeros(N_flybys, length_GAaF(flybyInd));
                vxTemp = zeros(N_flybys, length_GAaF(flybyInd));
                vyTemp = zeros(N_flybys, length_GAaF(flybyInd));
                vzTemp = zeros(N_flybys, length_GAaF(flybyInd));
                tTemp = zeros(N_flybys, length_GAaF(flybyInd));
                
                xTemp(1:flybyInd-1,1:maxGAaF) = GAa_xF(1:flybyInd-1,:);
                yTemp(1:flybyInd-1,1:maxGAaF) = GAa_yF(1:flybyInd-1,:);
                zTemp(1:flybyInd-1,1:maxGAaF) = GAa_zF(1:flybyInd-1,:);
                vxTemp(1:flybyInd-1,1:maxGAaF) = GAa_vxF(1:flybyInd-1,:);
                vyTemp(1:flybyInd-1,1:maxGAaF) = GAa_vyF(1:flybyInd-1,:);
                vzTemp(1:flybyInd-1,1:maxGAaF) = GAa_vzF(1:flybyInd-1,:);
                tTemp(1:flybyInd-1,1:maxGAaF) = GAa_tF(1:flybyInd-1,:);
                
                GAa_xF = xTemp;
                GAa_yF = yTemp;
                GAa_zF = zTemp;
                GAa_vxF = vxTemp;
                GAa_vyF = vyTemp;
                GAa_vzF = vzTemp;
                GAa_tF = tTemp;
            end
            
            % Store spacecraft state over flyby
            GAa_xF(flybyInd,1:length_GAaF(flybyInd)) = state(:,1);
            GAa_yF(flybyInd,1:length_GAaF(flybyInd)) = state(:,2);
            GAa_zF(flybyInd,1:length_GAaF(flybyInd)) = state(:,3);
            GAa_vxF(flybyInd,1:length_GAaF(flybyInd)) = state(:,4);
            GAa_vyF(flybyInd,1:length_GAaF(flybyInd)) = state(:,5);
            GAa_vzF(flybyInd,1:length_GAaF(flybyInd)) = state(:,6);
            GAa_tF(flybyInd,1:length_GAaF(flybyInd)) = ta;
            
            maxGAaF = max(length_GAaF);
            
         %% Backwards shooting, first GA phase
            
            X0 = [rp(:,flybyInd).', vp(:,flybyInd).'];

            tspan = [tEnd*TU/TU_flybys(bShort) + dt_GAa(flybyInd), tEnd*TU/TU_flybys(bShort) + dt_GAa(flybyInd)/2];
            
            [ta,state] = ode45(@(t,X0) EoM_hifi_coast(t, X0, mu_FBCent(flybyInd), N_perts, mu_pertFB(:,flybyInd), x_centPertBig, ...
                y_centPertBig, z_centPertBig, FB_impulse_timesa(i,:), TU_flybys(bShort)), tspan, X0);
            
            % Figure out how many data points are in phase now
            length_GAaB(flybyInd) = length(state(:,1));
            
            % Sort out any difference in size between different vectors
            if flybyInd == 1
                GAa_xB = zeros(N_flybys,length_GAaB(1));
                GAa_yB = zeros(N_flybys,length_GAaB(1));
                GAa_zB = zeros(N_flybys,length_GAaB(1));
                GAa_vxB = zeros(N_flybys,length_GAaB(1));
                GAa_vyB = zeros(N_flybys,length_GAaB(1));
                GAa_vzB = zeros(N_flybys,length_GAaB(1));
                GAa_tB = zeros(N_flybys,length_GAaB(1));
                
            % If more points in this state, need a bigger array for storage
            elseif length_GAaB(flybyInd) > maxGAaB
                xTemp = zeros(N_flybys, length_GAaB(flybyInd));
                yTemp = zeros(N_flybys, length_GAaB(flybyInd));
                zTemp = zeros(N_flybys, length_GAaB(flybyInd));
                vxTemp = zeros(N_flybys, length_GAaB(flybyInd));
                vyTemp = zeros(N_flybys, length_GAaB(flybyInd));
                vzTemp = zeros(N_flybys, length_GAaB(flybyInd));
                tTemp = zeros(N_flybys, length_GAaB(flybyInd));
                
                xTemp(1:flybyInd-1,1:maxGAaB) = GAa_xB(1:flybyInd-1,:);
                yTemp(1:flybyInd-1,1:maxGAaB) = GAa_yB(1:flybyInd-1,:);
                zTemp(1:flybyInd-1,1:maxGAaB) = GAa_zB(1:flybyInd-1,:);
                vxTemp(1:flybyInd-1,1:maxGAaB) = GAa_vxB(1:flybyInd-1,:);
                vyTemp(1:flybyInd-1,1:maxGAaB) = GAa_vyB(1:flybyInd-1,:);
                vzTemp(1:flybyInd-1,1:maxGAaB) = GAa_vzB(1:flybyInd-1,:);
                tTemp(1:flybyInd-1,1:maxGAaB) = GAa_tB(1:flybyInd-1,:);
                
                GAa_xB = xTemp;
                GAa_yB = yTemp;
                GAa_zB = zTemp;
                GAa_vxB = vxTemp;
                GAa_vyB = vyTemp;
                GAa_vzB = vzTemp;
                GAa_tB = tTemp;
            end
            
            % Store spacecraft state over flyby
            GAa_xB(flybyInd,1:length_GAaB(flybyInd)) = flipud(state(:,1));
            GAa_yB(flybyInd,1:length_GAaB(flybyInd)) = flipud(state(:,2));
            GAa_zB(flybyInd,1:length_GAaB(flybyInd)) = flipud(state(:,3));
            GAa_vxB(flybyInd,1:length_GAaB(flybyInd)) = flipud(state(:,4));
            GAa_vyB(flybyInd,1:length_GAaB(flybyInd)) = flipud(state(:,5));
            GAa_vzB(flybyInd,1:length_GAaB(flybyInd)) = flipud(state(:,6));
            GAa_tB(flybyInd,1:length_GAaB(flybyInd)) = flipud(ta);
            
            maxGAaB = max(length_GAaB);
            
            %% Second phase goes periapse back out to SOI
                        
            x_centPertBig = zeros(N_perts, N*2);
            y_centPertBig = zeros(N_perts, N*2);
            z_centPertBig = zeros(N_perts, N*2);
            for j = 1:N_perts
                x_centPertBig(j, :) = stateFBtoPertb(1+(j-1)*3+(flybyInd-1)*N_perts*3,:);
                y_centPertBig(j, :) = stateFBtoPertb(2+(j-1)*3+(flybyInd-1)*N_perts*3,:);
                z_centPertBig(j, :) = stateFBtoPertb(3+(j-1)*3+(flybyInd-1)*N_perts*3,:);
            end

            % Forwards shooting
            X0 = [rp(:,flybyInd).', vp(:,flybyInd).'];
%             tEnd*TU
%             dt_GAa(flybyInd)*TU_flybys(bShort)/86400
            tp = tEnd*TU/TU_flybys(bShort) + dt_GAa(flybyInd);
            tspan = [tp, tp + dt_GAb(flybyInd)/2];
            
            [tb,state] = ode45(@(t,X0) EoM_hifi_coast(t, X0, mu_FBCent(flybyInd), N_perts, mu_pertFB(:,flybyInd), x_centPertBig, ...
                y_centPertBig, z_centPertBig, FB_impulse_timesb(i,:), TU_flybys(bShort)), tspan, X0);
            
            % Figure out how many data points are in phase now
            length_GAbF(flybyInd) = length(state(:,1));
            
            % Sort out any difference in size between different vectors
            if flybyInd == 1
                GAb_xF = zeros(N_flybys,length_GAbF(1));
                GAb_yF = zeros(N_flybys,length_GAbF(1));
                GAb_zF = zeros(N_flybys,length_GAbF(1));
                GAb_vxF = zeros(N_flybys,length_GAbF(1));
                GAb_vyF = zeros(N_flybys,length_GAbF(1));
                GAb_vzF = zeros(N_flybys,length_GAbF(1));
                GAb_tF = zeros(N_flybys,length_GAbF(1));
                
            % If more points in this state, need a bigger array for storage
            elseif length_GAbF(flybyInd) > maxGAbF
                xTemp = zeros(N_flybys, length_GAbF(flybyInd));
                yTemp = zeros(N_flybys, length_GAbF(flybyInd));
                zTemp = zeros(N_flybys, length_GAbF(flybyInd));
                vxTemp = zeros(N_flybys, length_GAbF(flybyInd));
                vyTemp = zeros(N_flybys, length_GAbF(flybyInd));
                vzTemp = zeros(N_flybys, length_GAbF(flybyInd));
                tTemp = zeros(N_flybys, length_GAbF(flybyInd));
                
                xTemp(1:flybyInd-1,1:maxGAbF) = GAb_xF(1:flybyInd-1,:);
                yTemp(1:flybyInd-1,1:maxGAbF) = GAb_yF(1:flybyInd-1,:);
                zTemp(1:flybyInd-1,1:maxGAbF) = GAb_zF(1:flybyInd-1,:);
                vxTemp(1:flybyInd-1,1:maxGAbF) = GAb_vxF(1:flybyInd-1,:);
                vyTemp(1:flybyInd-1,1:maxGAbF) = GAb_vyF(1:flybyInd-1,:);
                vzTemp(1:flybyInd-1,1:maxGAbF) = GAb_vzF(1:flybyInd-1,:);
                tTemp(1:flybyInd-1,1:maxGAbF) = GAb_tF(1:flybyInd-1,:);
                
                GAb_xF = xTemp;
                GAb_yF = yTemp;
                GAb_zF = zTemp;
                GAb_vxF = vxTemp;
                GAb_vyF = vyTemp;
                GAb_vzF = vzTemp;
                GAb_tF = tTemp;
            end
            
            % Store spacecraft state over flyby
            GAb_xF(flybyInd,1:length_GAbF(flybyInd)) = state(:,1);
            GAb_yF(flybyInd,1:length_GAbF(flybyInd)) = state(:,2);
            GAb_zF(flybyInd,1:length_GAbF(flybyInd)) = state(:,3);
            GAb_vxF(flybyInd,1:length_GAbF(flybyInd)) = state(:,4);
            GAb_vyF(flybyInd,1:length_GAbF(flybyInd)) = state(:,5);
            GAb_vzF(flybyInd,1:length_GAbF(flybyInd)) = state(:,6);
            GAb_tF(flybyInd,1:length_GAbF(flybyInd)) = tb;
            
            maxGAbF = max(length_GAbF);
            
            %% Backwards shooting, second GA phase. From SOI backwards
            
            vEndGAb = vinff(:,flybyInd)*VU/VU_flybys(bShort);
            rEndGAb = radec2cart(SOI_flybys(bShort), RAf(flybyInd), DECf(flybyInd));
            
            % Forwards shooting
            X0 = [rEndGAb, vEndGAb.'];
            
            tspan = [tp + dt_GAb(flybyInd), tp + dt_GAb(flybyInd)/2];
            
            [tb,state] = ode45(@(t,X0) EoM_hifi_coast(t, X0, mu_FBCent(flybyInd), N_perts, mu_pertFB(:,flybyInd), x_centPertBig, ...
                y_centPertBig, z_centPertBig, FB_impulse_timesb(i,:), TU_flybys(bShort)), tspan, X0);
            
            % Figure out how many data points are in phase now
            length_GAbB(flybyInd) = length(state(:,1));
            
            % Sort out any difference in size between different vectors
            if flybyInd == 1
                GAb_xB = zeros(N_flybys,length_GAbB(1));
                GAb_yB = zeros(N_flybys,length_GAbB(1));
                GAb_zB = zeros(N_flybys,length_GAbB(1));
                GAb_vxB = zeros(N_flybys,length_GAbB(1));
                GAb_vyB = zeros(N_flybys,length_GAbB(1));
                GAb_vzB = zeros(N_flybys,length_GAbB(1));
                GAb_tB = zeros(N_flybys,length_GAbB(1));
                
            % If more points in this state, need a bigger array for storage
            elseif length_GAbB(flybyInd) > maxGAbB
                xTemp = zeros(N_flybys, length_GAbB(flybyInd));
                yTemp = zeros(N_flybys, length_GAbB(flybyInd));
                zTemp = zeros(N_flybys, length_GAbB(flybyInd));
                vxTemp = zeros(N_flybys, length_GAbB(flybyInd));
                vyTemp = zeros(N_flybys, length_GAbB(flybyInd));
                vzTemp = zeros(N_flybys, length_GAbB(flybyInd));
                tTemp = zeros(N_flybys, length_GAbB(flybyInd));
                
                xTemp(1:flybyInd-1,1:maxGAbB) = GAb_xB(1:flybyInd-1,:);
                yTemp(1:flybyInd-1,1:maxGAbB) = GAb_yB(1:flybyInd-1,:);
                zTemp(1:flybyInd-1,1:maxGAbB) = GAb_zB(1:flybyInd-1,:);
                vxTemp(1:flybyInd-1,1:maxGAbB) = GAb_vxB(1:flybyInd-1,:);
                vyTemp(1:flybyInd-1,1:maxGAbB) = GAb_vyB(1:flybyInd-1,:);
                vzTemp(1:flybyInd-1,1:maxGAbB) = GAb_vzB(1:flybyInd-1,:);
                tTemp(1:flybyInd-1,1:maxGAbB) = GAb_tB(1:flybyInd-1,:);
                
                GAb_xB = xTemp;
                GAb_yB = yTemp;
                GAb_zB = zTemp;
                GAb_vxB = vxTemp;
                GAb_vyB = vyTemp;
                GAb_vzB = vzTemp;
                GAb_tB = tTemp;
            end
            
            % Store spacecraft state over flyby
            GAb_xB(flybyInd,1:length_GAbB(flybyInd)) = flipud(state(:,1));
            GAb_yB(flybyInd,1:length_GAbB(flybyInd)) = flipud(state(:,2));
            GAb_zB(flybyInd,1:length_GAbB(flybyInd)) = flipud(state(:,3));
            GAb_vxB(flybyInd,1:length_GAbB(flybyInd)) = flipud(state(:,4));
            GAb_vyB(flybyInd,1:length_GAbB(flybyInd)) = flipud(state(:,5));
            GAb_vzB(flybyInd,1:length_GAbB(flybyInd)) = flipud(state(:,6));
            GAb_tB(flybyInd,1:length_GAbB(flybyInd)) = flipud(tb);
            
            maxGAbB = max(length_GAbB);
            
            %% Back to regular frame for next phase
            
            rEnd = rEndGAb*DU_flybys(bShort)/DU + [pert_FBxb(bCurrent, flybyInd*N*2), ...
                pert_FByb(bCurrent, flybyInd*N*2), pert_FBzb(bCurrent, flybyInd*N*2)]*DU_flybys(bShort)/DU;
            vEnd = flyby_vEnd(:,flybyInd) + vinff(:,flybyInd);
            vEnd = vEnd.';
            tEnd = (tp + dt_GAb(flybyInd))*TU_flybys(bShort)/TU;
        end
        flybyInd = flybyInd+1;
    end
    
%% Set up end control node as next phase start control node

    rStart = rEnd;
    vStart = vEnd;
    mStart = mEnd;
    tStart = tEnd;
    
end

%% Plot trajectory

figure(1)

plot(0,0,'ko')
hold on
j=1;
% Plot flyby bodies
for i = 1:NpCurrent
    if any(i == whichFlyby)
        flybyInd = flybySequence(j);
        plot3(pert_allx(flybyInd, i*N-(N-1):i*N)*DU, pert_ally(flybyInd, i*N-(N-1):i*N)*DU, pert_allz(flybyInd, i*N-(N-1):i*N)*DU, '-')
        j=j+1;
    end
end

if startBody > 2
    plot3(pert_allx(startBody-2, 1:N)*DU, pert_ally(startBody-2, 1:N)*DU, pert_allz(startBody-2, 1:N)*DU, '-')
end

if NpCurrent == Np && endBody > 2
    plot3(pert_allx(endBody-2, N*Np-(N-1):N*Np)*DU, pert_ally(endBody-2, N*Np-(N-1):N*Np)*DU, pert_allz(endBody-2, N*Np-(N-1):N*Np)*DU, '-')
end

for i = 1:NpCurrent
    plot3(x_all(i, 1:N/2+2)*DU, y_all(i, 1:N/2+2)*DU, z_all(i, 1:N/2+2)*DU, 'k-')
    plot3(x_all(i, N/2+3:end)*DU, y_all(i, N/2+3:end)*DU, z_all(i, N/2+3:end)*DU, 'g-')
end

xlabel('x direction [km]')
ylabel('y direction [km]')
zlabel('z direction [km]')
axis equal

%% Plot velocity

figure(2)
hold on

for i = 1:NpCurrent
    v = [vx_all(i,:); vy_all(i,:); vz_all(i,:)];
    plot(t_all(i,1:N/2+2)*TU, vecnorm(v(:,1:N/2+2), 2,1)*VU, 'kx')
    plot(t_all(i,N/2+3:end)*TU, vecnorm(v(:,N/2+3:end), 2,1)*VU, 'gx')
end

xlabel('Time since J2000 [s]')
ylabel('Velocity magnitude [km/s]')

%% Plot mass

figure(3)
hold on
thrustInd = 1;
for i = 1:NpCurrent
    if any(i == whichThrust)
        % Start mass is start node, end mass is impulse before match point
        plot(t_all(i,1:N/2+1)*TU, m_all(i,1:N/2+1)*MU, 'k-')
        % Start mass is impulse after match point, end is impulse before end
        plot(t_all(i,N/2+4:end)*TU, m_all(i,N/2+2:end)*MU, 'k-')
        thrustInd = thrustInd+1;
        lastMf = m_all(i,end);
    else
        if thrustInd > 1
            plot(t_all(i,:)*TU, lastMf*ones(1,N+4)*MU, 'k-')
        else
            plot(t_all(i,:)*TU, m0*ones(1,N+4)*MU, 'k-')
        end
    end
end

xlabel('Time since J2000 [s]')
ylabel('Mass [kg]')

%% Plot GAs

flybyInd = 1;
for i = 1:flybyPresent
    if FBEnd(i) == 1
        bCurrent = flybySequence(flybyInd);
        bShort = find(flybyBodies == bCurrent);
        figure(3+i)
        plot3(GAa_xF(flybyInd, 1:length_GAaF(flybyInd))*DU_flybys(bShort), GAa_yF(flybyInd, 1:length_GAaF(flybyInd))*DU_flybys(bShort), GAa_zF(flybyInd, 1:length_GAaF(flybyInd))*DU_flybys(bShort), 'k-')
        hold on
        plot3(GAa_xB(flybyInd, 1:length_GAaB(flybyInd))*DU_flybys(bShort), GAa_yB(flybyInd, 1:length_GAaB(flybyInd))*DU_flybys(bShort), GAa_zB(flybyInd, 1:length_GAaB(flybyInd))*DU_flybys(bShort), 'b-')
        plot3(GAb_xF(flybyInd, 1:length_GAbF(flybyInd))*DU_flybys(bShort), GAb_yF(flybyInd, 1:length_GAbF(flybyInd))*DU_flybys(bShort), GAb_zF(flybyInd, 1:length_GAbF(flybyInd))*DU_flybys(bShort), 'g-')
        plot3(GAb_xB(flybyInd, 1:length_GAbB(flybyInd))*DU_flybys(bShort), GAb_yB(flybyInd, 1:length_GAbB(flybyInd))*DU_flybys(bShort), GAb_zB(flybyInd, 1:length_GAbB(flybyInd))*DU_flybys(bShort), 'r-')
        plot3(0,0,0,'ko')
        plot3(rp(1,flybyInd)*DU_flybys(bShort), rp(2,flybyInd)*DU_flybys(bShort), rp(3,flybyInd)*DU_flybys(bShort), 'kx')
    end
    flybyInd = flybyInd + 1;
end

%% Plot GA velocity

flybyInd = 1;
for i = 1:flybyPresent
    if FBEnd(i) == 1
        bCurrent = flybySequence(flybyInd);
        bShort = find(flybyBodies == bCurrent);
        figure(3+flybyPresent+i)
        v = [GAa_vxF(flybyInd,1:length_GAaF(flybyInd)); GAa_vyF(flybyInd,1:length_GAaF(flybyInd)); GAa_vzF(flybyInd,1:length_GAaF(flybyInd))];
        plot(GAa_tF(flybyInd,1:length_GAaF(flybyInd))*TU_flybys(bShort), vecnorm(v(:,1:length_GAaF(flybyInd)), 2,1)*VU_flybys(bShort), 'k-')
        hold on
        v = [GAa_vxB(flybyInd,1:length_GAaB(flybyInd)); GAa_vyB(flybyInd,1:length_GAaB(flybyInd)); GAa_vzB(flybyInd,1:length_GAaB(flybyInd))];
        plot(GAa_tB(flybyInd,1:length_GAaB(flybyInd))*TU_flybys(bShort), vecnorm(v(:,1:length_GAaB(flybyInd)), 2,1)*VU_flybys(bShort), 'k-')
        v = [GAb_vxF(flybyInd,1:length_GAbF(flybyInd)); GAb_vyF(flybyInd,1:length_GAbF(flybyInd)); GAb_vzF(flybyInd,1:length_GAbF(flybyInd))];
        plot(GAb_tF(flybyInd,1:length_GAbF(flybyInd))*TU_flybys(bShort), vecnorm(v(:,1:length_GAbF(flybyInd)), 2,1)*VU_flybys(bShort), 'k-')
        v = [GAb_vxB(flybyInd,1:length_GAbB(flybyInd)); GAb_vyB(flybyInd,1:length_GAbB(flybyInd)); GAb_vzB(flybyInd,1:length_GAbB(flybyInd))];
        plot(GAb_tB(flybyInd,1:length_GAbB(flybyInd))*TU_flybys(bShort), vecnorm(v(:,1:length_GAbB(flybyInd)), 2,1)*VU_flybys(bShort), 'k-')
    end
    flybyInd = flybyInd + 1;
end

%% Constraints

noCeqs = 7*NpCurrent + 12*flybyPresent;
% If start or end on SOI, need extra constraint forcing start or endpoint
% to be on SOI
if startBody == 1
    noCeqs = noCeqs + 1;
end
if endBody == 1
    noCeqs = noCeqs + 1;
end
noCs = N*thrustPresent + 3*flybyPresent;
ceq = zeros(1,noCeqs);
c = zeros(1,noCs);

% Equality constraints: has form c = 0
% Match point constraints at every phase
for i = 1:NpCurrent
   ceq(1+(i-1)*6) = x_all(i, N/2+2) - x_all(i, N/2+3);
   ceq(2+(i-1)*6) = y_all(i, N/2+2) - y_all(i, N/2+3);
   ceq(3+(i-1)*6) = z_all(i, N/2+2) - z_all(i, N/2+3);
   ceq(4+(i-1)*6) = vx_all(i, N/2+2) - vx_all(i, N/2+3);
   ceq(5+(i-1)*6) = vy_all(i, N/2+2) - vy_all(i, N/2+3);
   ceq(6+(i-1)*6) = vz_all(i, N/2+2) - vz_all(i, N/2+3);
   ceqInd = 6+(i-1)*6+1;
end

cInd = 0;

% Mass equality constraints, and up-to-unit control vector, for thrust arcs
for i = 1:NpCurrent
    ceq(ceqInd) = m_all(i, N/2+1) - m_all(i, N/2+2);
    ceqInd = ceqInd+1;
end
for i = 1:thrustPresent
    u_this = u((1:3)+(i-1)*3,:);
    c(1+N*(i-1) : N*i) = vecnorm(u_this,2,1)-1;
    cInd = N*i;
end

% vinf constraint for each GA and periapse height check
for i = 1:flybyPresent
    if FBEnd(i) == 1
        bShort = find(flybyBodies == flybySequence(i));
        ceq(ceqInd) = GAa_xF(i, end) - GAa_xB(i,1);
        ceq(ceqInd+1) = GAa_yF(i, end) - GAa_yB(i,1);
        ceq(ceqInd+2) = GAa_zF(i, end) - GAa_zB(i,1);
        ceq(ceqInd+3) = GAa_vxF(i, end) - GAa_vxB(i,1);
        ceq(ceqInd+4) = GAa_vyF(i, end) - GAa_vyB(i,1);
        ceq(ceqInd+5) = GAa_vzF(i, end) - GAa_vzB(i,1);
        
        ceq(ceqInd+6) = GAb_xF(i, end) - GAb_xB(i,1);
        ceq(ceqInd+7) = GAb_yF(i, end) - GAb_yB(i,1);
        ceq(ceqInd+8) = GAb_zF(i, end) - GAb_zB(i,1);
        ceq(ceqInd+9) = GAb_vxF(i, end) - GAb_vxB(i,1);
        ceq(ceqInd+10) = GAb_vyF(i, end) - GAb_vyB(i,1);
        ceq(ceqInd+11) = GAb_vzF(i, end) - GAb_vzB(i,1);
        ceqInd = ceqInd + 12;
        
        c(cInd) = h_mins(i)*DU/DU_flybys(bShort) + r_flybys(i)*DU/DU_flybys(bShort) - norm(rp(:,i));
        c(cInd+1) = norm(rp(:,i)) - h_maxs(i)*DU/DU_flybys(bShort) - r_flybys(i)*DU/DU_flybys(bShort);
        cInd = cInd+2;
    else
        ceq(ceqInd) = norm(vinfi(:,i)) - norm(vinff(:,i));
        c(cInd) = h_mins(i) + r_flybys(i) - rp_all(i);
        c(cInd+1) = rp_all(i) - h_maxs(i) - r_flybys(i);
        ceqInd = ceqInd+2;
    end
    
end

if startBody == 1
    ceq(ceqInd) = norm(ri) - SOI;
    ceqInd = ceqInd + 1;
end

if NpCurrent == Np && endBody == 1
    ceq(ceqInd) = norm(rf) - SOI;
end
end
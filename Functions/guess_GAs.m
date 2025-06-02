function guess = guess_GAs(x, consts)

% Constants
Np = consts(1);         % Number of phases
N = consts(2);          % Number of segments per phase
N_perts = consts(3);            % Number of perturbing bodies
N_ephem = consts(4);            % Number of point in ephemeris data
N_thrust = consts(5);           % Number of phases with thrust in
N_flybys = consts(6);
NpCurrent = consts(7);          % What number phase currently in
mu_central = consts(8);         % Central body standard gravitational parameter in DU^3/TU^2
% r_central = consts(9);          % Central body radius in DU
SOI = consts(10);                % Central body sphere of influence radius in DU
% m0 = consts(11);                % Initial spacecraft mass in MU
% N_act = consts(12);             % Number of active thrusters
% D = consts(13);                 % Thrust duty cycle
% Tmax = consts(14);              % Max thrust in MU DU/TU^2 (from kN)
% Isp = consts(15);               % Specific impulse in TU
% g0 = consts(16);                % Acceleration due to gravity on Earth's surface in DU/TU^2
DU = consts(17);                % Distance unit in km
TU = consts(18);                % Time unit in s
% MU = consts(19);                % Mass unit in kg
whichThrust = consts(20 : 19+N_thrust);     % Index of phases with thrust
flybySequence = consts(20+N_thrust : 19+N_thrust+N_flybys);   % Sequence of flybys
startBody = consts(20+N_thrust+N_flybys);
endBody = consts(21+N_thrust+N_flybys);
whichFlyby = consts(22+N_thrust+N_flybys : 21+N_thrust+2*N_flybys);
mu_perts = consts(22+N_thrust+2*N_flybys : 21+N_thrust+2*N_flybys+N_perts); % Standard gravitational parameter of perturbing bodies in DU^3/TU^2
% r_flybys = consts(22+N_thrust+2*N_flybys+N_perts : 21+N_thrust+3*N_flybys+N_perts);   % Radius of flyby bodies in DU
% h_mins = consts(22+N_thrust+3*N_flybys+N_perts : 21+N_thrust+4*N_flybys+N_perts);   % Min flyby heights in DU
% h_maxs = consts(22+N_thrust+4*N_flybys+N_perts : 21+N_thrust+5*N_flybys+N_perts);   % Max flyby heights in DU
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
% MU_flybys = consts(23+N_thrust+5*N_flybys+N_perts+(7+6*(N_perts-1))*N_ephem+3*N_fbodies : 22+N_thrust+5*N_flybys+N_perts+(7+6*(N_perts-1))*N_ephem+4*N_fbodies);

VU = DU/TU;
VU_flybys = DU_flybys./TU_flybys;
% thrustPresent = length(find(whichThrust<=NpCurrent,1,'last'));
flybyPresent = find(whichFlyby<=NpCurrent,1,'last');

%% Decision variables

% For initial phase: initial epoch, initial RA DEC, initial vinf relative to Jupiter
% For each phase: tof, mf, vinfi and vinff for the gravity assist at the
% end of the trajectory, control vector u
% For final phase: vinfi but no vinff because not doing a GA

dt = zeros(1, Np);
vinfi = zeros(3, N_flybys);
vinff = zeros(3, N_flybys);
mf = zeros(1, N_thrust);
u = zeros(3*N_thrust, N);

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
        vinfi(:,flybyInd) = x(decisionInd+1 : decisionInd+3);
        vinff(:,flybyInd) = x(decisionInd+4 : decisionInd+6);
        decisionInd = decisionInd+7;
        flybyInd = flybyInd + 1;
        
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
end

% freePresent = freeInd - 1;

%% Spline for ephemeris data

% Need a position for each perturbing body at the centre of each segment.
% Need velocity of each flyby body at the control node where the GA is
impulse_times = zeros(NpCurrent, N);
pert_allx = zeros(N_perts, NpCurrent*N);
pert_ally = zeros(N_perts, NpCurrent*N);
pert_allz = zeros(N_perts, NpCurrent*N);
flyby_v = zeros(3, N_flybys);
tstart = t0*TU;

for i = 1:NpCurrent
    thalf = (dt(i)/(2*N))*TU;
    tend = tstart + dt(i)*TU;
    impulse_times(i,:) = linspace(tstart+thalf, tend-thalf, N);
    tstart = tstart + dt(i)*TU;
    
    % Position values for each perturbing body at each impulse
    for j = 1:N_perts
        pert_allx(j, i*N-(N-1):i*N) = spline(et, stateLong(1+(j-1)*6, :), impulse_times(i,:));
        pert_ally(j, i*N-(N-1):i*N) = spline(et, stateLong(2+(j-1)*6, :), impulse_times(i,:));
        pert_allz(j, i*N-(N-1):i*N) = spline(et, stateLong(3+(j-1)*6, :), impulse_times(i,:));
    end
end

% Get velocities for the flybys
for i = 1:flybyPresent
    index = flybySequence(i);
    epoch = (t0 + sum(dt(1:i)))*TU;
    flyby_v(:, i) = spline(et, stateLong((4:6)+(index-1)*6,:), epoch);
end

%% Find periapse state for each flyby

rpMag = zeros(1,N_flybys);
rp = zeros(3, N_flybys);
vpMag = zeros(1, N_flybys);
vp = zeros(3, N_flybys);
dt_GAa = zeros(1, N_flybys);
dt_GAb = zeros(1, N_flybys);
xGAa = zeros(1, N_flybys);
yGAa = zeros(1, N_flybys);
zGAa = zeros(1, N_flybys);
vxGAa = zeros(1, N_flybys);
vyGAa = zeros(1, N_flybys);
vzGAa = zeros(1, N_flybys);
xGAb = zeros(1, N_flybys);
yGAb = zeros(1, N_flybys);
zGAb = zeros(1, N_flybys);
vxGAb = zeros(1, N_flybys);
vyGAb = zeros(1, N_flybys);
vzGAb = zeros(1, N_flybys);
RAi = zeros(1, N_flybys);
DECi = zeros(1, N_flybys);
RAf = zeros(1, N_flybys);
DECf = zeros(1, N_flybys);
vinfiProp = zeros(3, N_flybys);
vinffProp = zeros(3, N_flybys);

for i = 1:N_flybys
    p = whichFlyby(i);          % Phase index of the flyby
    bCurrent = flybySequence(i);
    bShort = find(flybyBodies == bCurrent);
    
    delta = acos(dot(vinfi(:,i), vinff(:,i))/(norm(vinfi(:,i))*norm(vinff(:,i))));
    e = 1/sin(delta/2);
    rphat = (vinfi(:,i) - vinff(:,i))/norm(vinfi(:,i) - vinff(:,i));       % might be abs
    hhat = cross(vinfi(:,i), vinff(:,i))/norm(cross(vinfi(:,i), vinff(:,i)));  % might be abs
    rpMag(i) = (mu_perts(bCurrent)*(e - 1))/(norm(vinfi(:,i))*norm(vinff(:,i)));       % dot used to be norm(vinf)
    rp(:,i) = rpMag(i)*rphat;

    vphat = cross(hhat, rphat);
    vpMag(i) = sqrt((norm(vinfi(:,i))*norm(vinff(:,i))) + (2*mu_perts(bCurrent))/rpMag(i));         % dot was norm
    vp(:,i) = vpMag(i)*vphat;
    
    %% Convert to the flyby-centred frame
    
    rpMag(i) = rpMag(i)*DU/DU_flybys(bShort);
    rp(:,i) = rp(:,i)*DU/DU_flybys(bShort);
    vpMag(i) = vpMag(i)*VU/VU_flybys(bShort);
    vp(:,i) = vp(:,i)*VU/VU_flybys(bShort);
    
    stateLong = stateLong*DU/DU_flybys(bShort);
    % SOI already in flyby centred, as dictated what DU is
    mu_central = mu_central*(DU^3/TU^2)/(DU_flybys(bShort)^3/TU_flybys(bShort)^2);
    mu_perts = mu_perts.*(DU^3/TU^2)./(DU_flybys(bShort)^3/TU_flybys(bShort)^2);
    
    stateFBcent = zeros(3,N_ephem);
    stateFBtoPert = zeros(3*N_perts,N_ephem);
    mu_pertsNew = zeros(1, N_perts);
    
    stateFBcent(1,:) = -stateLong(1+(bCurrent-1)*6,:);
    stateFBcent(2,:) = -stateLong(2+(bCurrent-1)*6,:);
    stateFBcent(3,:) = -stateLong(3+(bCurrent-1)*6,:);
    
    for j = 1:N_perts        
        if j ~= bCurrent
            % These are all relative to the previous central body. To get them
            % relative to flyby body:
            % Vector from flyby body to central body, as central body is now
            % perturbation

            % Vector from flyby body to the other perturbing bodies
            stateFBtoPert(1+(j-1)*3,:) = stateFBcent(1,:) + stateLong(1+(j-1)*6,:);
            stateFBtoPert(2+(j-1)*3,:) = stateFBcent(2,:) + stateLong(2+(j-1)*6,:);
            stateFBtoPert(3+(j-1)*3,:) = stateFBcent(3,:) + stateLong(3+(j-1)*6,:);
            mu_pertsNew(j) = mu_perts(j);
        end
    end
    
    % Remove the empty gap where the new central body was in perturbation
    % state
    stateFBtoPert = nonzeros(stateFBtoPert);
    stateFBtoPert = reshape(stateFBtoPert, 3*(N_perts-1),N_ephem);  
    mu_pertsNew = nonzeros(mu_pertsNew).';      % Need to flip because nonzeros makes it a column vector
    
    mu_centNew = mu_perts(bCurrent);
    mu_pertsNew = [mu_pertsNew, mu_central];
    
    epoch = (t0 + sum(dt(1:p)))*TU/TU_flybys(bShort);
    
    %% Use a root finder to get the tof of each GA
    
    stateR = zeros(1,N_perts*N_ephem*3);
    for j = 1:N_perts-1
        stateR((1+3*(j-1)-1)*N_ephem + 1 : (1+3*(j-1))*N_ephem) = stateFBtoPert(1+(j-1)*3,:);
        stateR((2+3*(j-1)-1)*N_ephem + 1 : (2+3*(j-1))*N_ephem) = stateFBtoPert(2+(j-1)*3,:);
        stateR((3+3*(j-1)-1)*N_ephem + 1 : (3+3*(j-1))*N_ephem) = stateFBtoPert(3+(j-1)*3,:);
    end
    
    % Last state is for the old central body
    stateR((1+3*(N_perts-1)-1)*N_ephem + 1 : (1+3*(N_perts-1))*N_ephem) = stateFBcent(1,:);
    stateR((2+3*(N_perts-1)-1)*N_ephem + 1 : (2+3*(N_perts-1))*N_ephem) = stateFBcent(2,:);
    stateR((3+3*(N_perts-1)-1)*N_ephem + 1 : (3+3*(N_perts-1))*N_ephem) = stateFBcent(3,:);
    
    root_consts = [N_perts, N_ephem, N_int, SOI_flybys(bShort), rp(:,i).', vp(:,i).', epoch, TU_flybys(bShort), mu_centNew, ...
        mu_pertsNew, et, stateR];
    
    dt_guess = (SOI_flybys(bShort)/vpMag(i));
    dt_GAa(i) = abs(fzero(@(x)does_tof_reach_distance(x, root_consts), dt_guess));
    dt_guess = -(SOI_flybys(bShort)/vpMag(i));
    dt_GAb(i) = abs(fzero(@(x)does_tof_reach_distance(x, root_consts), dt_guess));
    
    %% Propagate backwards from periapse to incoming SOI point
    
    x_centPertBig = zeros(N_perts, N_int);
    y_centPertBig = zeros(N_perts, N_int);
    z_centPertBig = zeros(N_perts, N_int);
    times = linspace(epoch-dt_GAa(i)*2, epoch+dt_GAb(i)*2,N_int)*TU_flybys(bShort);
    for j = 1:N_perts
        x_centPertBig(j, :) = spline(et, stateR((1+3*(j-1)-1)*N_ephem + 1 : (1+3*(j-1))*N_ephem), times);
        y_centPertBig(j, :) = spline(et, stateR((2+3*(j-1)-1)*N_ephem + 1 : (2+3*(j-1))*N_ephem), times);
        z_centPertBig(j, :) = spline(et, stateR((3+3*(j-1)-1)*N_ephem + 1 : (3+3*(j-1))*N_ephem), times);
    end
    
    X0 = [rp(:,i).', vp(:,i).'];
    
    tspan = [epoch, epoch - dt_GAa(i)];
    
    [ta,state] = ode45(@(t,X0) EoM_hifi_coast(t, X0, mu_centNew, N_perts, mu_pertsNew, x_centPertBig, ...
        y_centPertBig, z_centPertBig, times, TU_flybys(bShort)), tspan, X0);
    
    xGAa(i) = state(end,1);
    yGAa(i) = state(end,2);
    zGAa(i) = state(end,3);
    vxGAa(i) = state(end,4);
    vyGAa(i) = state(end,5);
    vzGAa(i) = state(end,6);
    
    %% Plotting
    
%     [~, RAi(i), DECi(i)] = cart2radec([xGAa(i),yGAa(i),zGAa(i)]);
%     figure(i)
%     hold on
%     plot3(state(:,1)*DU_flybys(bShort), state(:,2)*DU_flybys(bShort), state(:,3)*DU_flybys(bShort), 'k--')
%     
%     % Test start position
%     ri = radec2cart(SOI_flybys(bShort), RAi(i), DECi(i));
%     
%     plot3(ri(1)*DU_flybys(bShort), ri(2)*DU_flybys(bShort), ri(3)*DU_flybys(bShort), 'gx')
% %     plot3(state(end,1)*DU_flybys(bShort), state(end,2)*DU_flybys(bShort), state(end,3)*DU_flybys(bShort), 'rx')
%     
%     figure(i+N_flybys)
%     hold on
%     vGA = [state(:,4), state(:,5), state(:,6)];
%     plot(ta*TU_flybys(bShort), vecnorm(vGA,2,2)*VU_flybys(bShort), 'k-')
    
    %% Propagate forwards from periapse to outgoing SOI point
    
    tspan = [epoch, epoch + dt_GAb(i)];
    
    [tb, state] = ode45(@(t,X0) EoM_hifi_coast(t, X0, mu_centNew, N_perts, mu_pertsNew, ...
        x_centPertBig, y_centPertBig, z_centPertBig, times, TU_flybys(bShort)), tspan, X0);
    
    xGAb(i) = state(end,1);
    yGAb(i) = state(end,2);
    zGAb(i) = state(end,3);
    vxGAb(i) = state(end,4);
    vyGAb(i) = state(end,5);
    vzGAb(i) = state(end,6);
    
    %% Plot again
    
%     figure(i)
%     hold on
%     [~, RAf(i), DECf(i)] = cart2radec([xGAb(i),yGAb(i),zGAb(i)]);
%     plot3(state(:,1)*DU_flybys(bShort), state(:,2)*DU_flybys(bShort), state(:,3)*DU_flybys(bShort), 'k--')
% %     plot3(rp(1,i)*DU_flybys(bShort), rp(2,i)*DU_flybys(bShort), rp(3,i)*DU_flybys(bShort), 'rx')
%     
%     % Test start position
%     rf = radec2cart(SOI_flybys(bShort), RAf(i), DECf(i));
%     plot3(rf(1)*DU_flybys(bShort), rf(2)*DU_flybys(bShort), rf(3)*DU_flybys(bShort), 'gx')
%     axis equal
% %     plot3(state(end,1)*DU_flybys(bShort), state(end,2)*DU_flybys(bShort), state(end,3)*DU_flybys(bShort), 'rx')
%     
%     th = 0:pi/50:2*pi;
%     xunit = SOI_flybys(bShort)*DU_flybys(bShort)* cos(th);
%     yunit = SOI_flybys(bShort)*DU_flybys(bShort)* sin(th);
%     plot(xunit, yunit, 'g--');
%     
%     figure(i+N_flybys)
%     hold on
%     vGA = [state(:,4), state(:,5), state(:,6)];
%     plot(tb*TU_flybys(bShort), vecnorm(vGA,2,2)*VU_flybys(bShort), 'k-')
    
    %% Convert to previous central body-centred frame
    
    stateLong = stateLong*DU_flybys(bShort)/DU;
%     SOI_flybys = SOI_flybys*DU_flybys(bShort)/DU;
    mu_central = mu_central*(DU_flybys(bShort)^3/TU_flybys(bShort)^2)/(DU^3/TU^2);
    mu_perts = mu_perts.*(DU_flybys(bShort)^3/TU_flybys(bShort)^2)./(DU^3/TU^2);
    
    % Points on SOI are in RA DEC
    [~, RAi(i), DECi(i)] = cart2radec([xGAa(i),yGAa(i),zGAa(i)]);
    [~, RAf(i), DECf(i)] = cart2radec([xGAb(i),yGAb(i),zGAb(i)]);
    
    % Velocities on SOI become vinfs of phases either side, so need in VU
    vinfiProp(:,i) = [vxGAa(i), vyGAa(i), vzGAa(i)]*VU_flybys(bShort)/VU;
    vinffProp(:,i) = [vxGAb(i), vyGAb(i), vzGAb(i)]*VU_flybys(bShort)/VU;
    
%     if i ==1
%         figure(5)
%         subplot(3,2,1)
%         plot(times, x_centPertBig(1,:), 'k-')
%         hold on
% 
%         subplot(3,2,2)
%         plot(times, y_centPertBig(1,:), 'k-')
%         hold on
% 
%         subplot(3,2,3)
%         plot(times, z_centPertBig(1,:), 'k-')
%         hold on
%         
%         subplot(3,2,4)
%         plot(times, x_centPertBig(2,:), 'k-')
%         hold on
% 
%         subplot(3,2,5)
%         plot(times, y_centPertBig(2,:), 'k-')
%         hold on
% 
%         subplot(3,2,6)
%         plot(times, z_centPertBig(2,:), 'k-')
%         hold on
%     end
end

%% Put stuff into guess

for i = 1:Np
    % If flyby at the end, deduct flyby tof from phase tof
    if any(i == whichFlyby)
        bShort = find(flybyBodies == flybySequence(i));
        dt(i) = dt(i) - dt_GAa(i)*TU_flybys(bShort)/TU;
    end
    
    % If had a flyby at the start, deduct flyby tof
    if any((i-1) == whichFlyby)
        bShort = find(flybyBodies == flybySequence(i-1));
        dt(i) = dt(i) - dt_GAb(i-1)*TU_flybys(bShort)/TU;
    end
end

guess = [RAi, DECi, RAf, DECf, vinfiProp(1,:),vinfiProp(2,:),vinfiProp(3,:), vinffProp(1,:),vinffProp(2,:),vinffProp(3,:), ...
    rp(1,:),rp(2,:),rp(3,:), vp(1,:),vp(2,:),vp(3,:), dt_GAa, dt_GAb, dt];

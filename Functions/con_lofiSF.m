function [c, ceq] = con_lofiSF(x,consts)

% Constants
Np = consts(1);                 % Number of phases
N = consts(2);                  % Number of segments per phase
N_perts = consts(3);            % Number of perturbing bodies
N_ephem = consts(4);            % Number of point in ephemeris data
N_thrust = consts(5);           % Number of phases with thrust in
N_flybys = consts(6);           % Number of GAs
NpCurrent = consts(7);          % What number phase currently in
mu_central = consts(8);         % Central body standard gravitational parameter in DU^3/TU^2
r_central = consts(9);          % Central body radius in DU
SOI = consts(10);               % Central body sphere of influence radius in DU
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
% If there are no flybys, still pass in a zero value for sequence and which
% flyby, so need to include a 1 value for the indexing
if N_flybys == 0
    FB_ind = 1;
else
    FB_ind = N_flybys;
end
flybySequence = consts(20+N_thrust : 19+N_thrust+FB_ind);   % Sequence of flybys
startBody = consts(20+N_thrust+FB_ind);                     % Where the trajectory starts
endBody = consts(21+N_thrust+FB_ind);                       % Where the trajectory ends
constInd = 22+N_thrust+FB_ind;
whichFlyby = consts(constInd : constInd+FB_ind-1);    % Which phases end in gravity assists
constInd = constInd+FB_ind;
mu_perts = consts(constInd : constInd+N_perts-1);       % Standard gravitational parameter of perturbing bodies in DU^3/TU^2
constInd = constInd + N_perts;
r_flybys = consts(constInd : constInd+FB_ind-1);      % Radius of flyby bodies in DU
constInd = constInd + FB_ind;
h_mins = consts(constInd : constInd+FB_ind-1);        % Min flyby heights in DU
constInd = constInd + FB_ind;
h_maxs = consts(constInd : constInd+FB_ind-1);        % Max flyby heights in DU
constInd = constInd + FB_ind;
et = consts(constInd : constInd+N_ephem-1);   % Epochs at which ephemerides taken in DU

% Extract ephemeris data
stateLong = zeros(6*N_perts, N_ephem);
for i = 1:N_perts
    for j = 1:6
        stateLong(j+(i-1)*6,:) = consts(constInd+(j+6*(i-1))*N_ephem : constInd-1+(j+6*(i-1)+1)*N_ephem);
    end
end

VU = DU/TU;                             % Velocity unit in km/s
thrustPresent = length(find(whichThrust<=NpCurrent));  % How many phases at this point have thrusting
flybyPresent = find(whichFlyby<=NpCurrent,1,'last');            % How many phases at this point end in a gravity assist

%% Decision variables

% For initial phase: initial epoch, initial RA DEC, initial vinf relative to Jupiter
% For each phase: tof, mf, vinfi and vinff for the gravity assist at the
% end of the trajectory, control vector u
% For final phase: vinfi but no vinff because not doing a GA

% Some empty arrays to put stuff in
dt = zeros(1, Np);
vinfi = zeros(3, N_flybys);
vinff = zeros(3, N_flybys);
mf = zeros(1, N_thrust);
u = zeros(3*N_thrust, N);

t0 = x(1);                                      % Initial epoch in TU
% If trajectory starts on SOI
if startBody == 1
    RA = x(2);                                  % Initial RA in rad
    DEC = x(3);                                 % Initial DEC in rad
    vi = x(4:6);                                % Initial velocity in VU
    ri = radec2cart(SOI, RA, DEC);              % Convert RA and DEC to point on SOI in Cartesian coordinates
    decisionInd = 7;                            % This keeps tabs on the problem size
    
    % Or, if trajectory starts in orbit around central body
elseif startBody == 2
    ai = x(2);                                  % Semimajor axis in DU
    AoPi = x(3);                                % Argument of periapse in rad
    ei = x(4);                                  % Eccentricity
    inci = x(5);                                % Inclination in rad
    RAANi = x(6);                               % Right ascension of ascending node in rad
    TAi = x(7);                                 % True anomaly in rad
    orbels = [ai, AoPi, ei, inci, RAANi, TAi];
    [ri, vi] = orb2state(mu_central, orbels);   % Convert orbital elements to position and velocity in Cartesian coordinates
    decisionInd = 8;
    
    % Or, if trajectory starts at location of perturbing body
else
    % Initial position and velocity in DU and VU are the position and
    % velocity of the body started at
    ri = spline(et, stateLong((1:3)+((startBody-2)-1)*6, :), t0*TU);
    vi = spline(et, stateLong((4:6)+((startBody-2)-1)*6, :), t0*TU).' + x(2:4);
    decisionInd = 5;
end

% Keep track of the number of thrust arcs, free points, and gravity assists
thrustInd = 1;
freeInd = 1;
flybyInd = 1;

% If there are any non-flyby intermediate control nodes, made some position and
% velocity arrays for them
if N_flybys ~= Np-1
    r_free = zeros(3,Np-1-N_flybys);
    v_free = zeros(3,Np-1-N_flybys);
end

% Go through the decision variables every phase has, and the decision
% variables for the last phase
for i = 1:NpCurrent
    dt(i) = x(decisionInd);                     % Time of flight of this phase in TU
    
    % If phase ends in a flyby
    if any(i == whichFlyby)
        % Need the initial and final gravity assist velocities relative to
        % flyby body, in VU
        vinfi(:,flybyInd) = x(decisionInd+1 : decisionInd+3);
        vinff(:,flybyInd) = x(decisionInd+4 : decisionInd+6);
        decisionInd = decisionInd+7;
        flybyInd = flybyInd + 1;
        
        % If it's the final phase
    elseif i == Np
        % If end at SOI
        if endBody == 1
            RAend = x(decisionInd+1);            % Final RA in rad
            DECend = x(decisionInd+2);           % Final DEC in rad
            vf = x(decisionInd+3:decisionInd+5); % Final velocity in VU
            rf = radec2cart(SOI, RAend, DECend); % Convert RA and DEC to position in DU with magnitude SOI
            decisionInd = decisionInd+6;    
            
            % Or, if trajectory ends in orbit around central body
        elseif endBody == 2
            af = x(decisionInd+1);               % Final semimajor axis in DU
            AoPf = x(decisionInd+2);             % Final argument of periapse in rad
            ef = x(decisionInd+3);               % Final eccentricity
            incf = x(decisionInd+4);             % Final inclination in rad
            RAANf = x(decisionInd+5);            % Final right ascension of ascending node in rad
            TAf = x(decisionInd+6);              % Final true anomaly in rad
            orbels = [af, AoPf, ef, incf, RAANf, TAf];
            [rf, vf] = orb2state(mu_central, orbels);   % Convert to final position and velocity in DU and VU
            decisionInd = decisionInd+7;

            % Or, if trajectory ends at location of perturbing body
        else
            % Get final position and velocity in DU and VU from the
            % position and velocity of the target body
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
    tstart = tend;      % Start time for next phase is end time of this phase
    
    % Position values for each perturbing body at each impulse
    for j = 1:N_perts
        pert_allx(j, i*N-(N-1):i*N) = spline(et, stateLong(1+(j-1)*6, :), impulse_times(i,:));
        pert_ally(j, i*N-(N-1):i*N) = spline(et, stateLong(2+(j-1)*6, :), impulse_times(i,:));
        pert_allz(j, i*N-(N-1):i*N) = spline(et, stateLong(3+(j-1)*6, :), impulse_times(i,:));
    end
end

% Get velocities for the flybys
if N_flybys > 0
    for i = 1:flybyPresent
        index = flybySequence(i);
        epoch = (t0 + sum(dt(1:i)))*TU;
        flyby_v(:, i) = spline(et, stateLong((4:6)+(index-1)*6,:), epoch);
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
        rEnd = [pert_allx(flybySequence(flybyInd), i*N), pert_ally(flybySequence(flybyInd), i*N), pert_allz(flybySequence(flybyInd), i*N)];
        vEnd = flyby_v(:,flybyInd) + vinfi(:,flybyInd);
        
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
        flybyInd = flybyInd+1;
    end
    
%% Set up end control node as next phase start control node

    rStart = rEnd;
    vStart = vEnd;
    mStart = mEnd;
    tStart = tEnd;
    
end

%% Constraints


noCeqs = 6*NpCurrent + thrustPresent + flybyPresent;
% If start or end on SOI, need extra constraint forcing start or endpoint
% to be on SOI
if startBody == 1
    noCeqs = noCeqs + 1;
end
if endBody == 1
    noCeqs = noCeqs + 1;
end
noCs = N*thrustPresent + flybyPresent*2;
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
    cInd = N*i+1;
end

% vinf constraint for each GA and periapse height check
if N_flybys > 0
    for i = 1:flybyPresent
        ceq(ceqInd) = norm(vinfi(:,i)) - norm(vinff(:,i));
        ceqInd = ceqInd+1;
        
        c(cInd) = h_mins(i) + r_flybys(i) - rp_all(i);
        c(cInd+1) =  rp_all(i) + r_flybys(i)  - h_maxs(i);
        cInd = cInd + 2;
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
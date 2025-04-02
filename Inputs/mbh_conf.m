
% Control vector. Need one in the x,y,z directions for each impulse in each
% phase that permits thrust. Shouldn't need to use anything besides bounds
% of -1 and 1.
u_min = (-1)*ones(N_thrust*3, N);
u_max = ones(N_thrust*3, N);

%% User Inputs - MBH parameters
% Not usually much need to change the defaults
% ^^THIS IS NOT TRUE, MBH PARAMETERS MUST BE CHANGED FOR DIFFERENT SYSTEMS

use_mbh = true
MBH_noLoops > 0 

t0Hop = 9*86400/TU;                % Max amount to hop launch epoch in TU
dtHop = 0.2*86400/TU;                 % Max amount to hop time of flights in TU
MBH_noLoops = 100;                  % Number of times to run
rho_hop = 0.2;                      % Probability of a hop
MBH_tail = 0.8;                     % MBH tail parameter
MBH_theta = 1;                      % MBH threshold/location parameter

% MBH scale parameters sigma
s_angles = 0.1;                     % Angles in rad
s_a = 100/DU;                       % Semimajor axis
s_e = 0.01;                         % Eccentricity
s_t0 = (10*86400)/TU;               % Initial epoch
s_dt = (0.1*86400)/TU;              % Time of flight
s_r = 100/DU;                       % Positions
s_v = 0.1/VU;                       % Velocities
s_u = 0.1;                          % Control vector
s_m = 10/MU;                        % Mass
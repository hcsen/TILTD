%% User inputs - mission constants

% Celestial system
N_perts = 2;                         % Number of perturbing bodies
centralBodyName = 'Saturn';             % Central body name in words
pertNames = strsplit('Enceladus, Dione', ',');       % Perturbing bodies in words
pertRelativeTo = 'SATURN BARYCENTER';     % System barycentre name in words

Np = 2;                              % Number of trajectory phases
N_flybys = 1;                        % Number of flybys

% Flybys numbered 1,2,3 etc according to order given in pertNames
flybySequence = [1];
whichFlyby = [1];                    % List which phases end in a flyby

% For start and end, 1 = SOI, 2 = orbit around central body, then index of the
% perturbing bodies + 2. Start at SOI specified with RA and DEC, orbit
% specified by Keplerian orbital elements a, AoP, e, inc, RAAN, t_anom.
startBody = 4;
endBody = 4;

N = 50;                              % Number of segments per trajectory phase

N_ephem = 2000;                      % Number of points over which to get ephemeris data
N_thrust = 2;                        % Number of thrust arcs
whichThrust = [1,2];                 % List of which phases have thrust

% Spacecraft info
N_act = 1;                           % Number of active thrusters
D = 1;                               % Thruster duty cycle = percentage of phase turned on
Tmax = 1.1e-6;                        % Max thrust in kN
Isp = 1700;                          % Isp in s
m0 = 14;                            % Initial mass in kg

% Tolerances for optimisation
NLP_feas_tol = 1e-4;                 % NLP solver feasibility tolerance
NLP_tol = 1e-1;                      % NLP solver optimality tolerance
NLP_steptol = 1e-10;                 % Max NLP solver runtime
NLP_iter_max = 10000;                % Max number of NLP solver major iterations

%% User inputs - Constants

G = 6.674e-20;                       % Gravitational constant in km^3 kg^-1 s^-2

% Provide masses of central body all perturbing bodies in order given
% before in kg
m_central = 5.683e26;                % Central body mass (Saturn)
m_perts = [1.080e20, 1.095e21];      % Enceladus, Dione

% Compute standard gravitational parameter of all perturbing bodies
mu_central = m_central*G;            % Central body standard gravitational parameter in km^3/s^2
mu_perts = m_perts.*G;               % Perturbing bodies standard gravitational parameter in km^3/s^2

% Radii of central body and perturbing bodies in km
r_central = 58232;                  % Radius of central body in km
r_perts = [252.1, 561.4];            % Enceladus, Dione

% Get sphere of influence radii
d_Scent = 1.434e9;                  % Average distance of Sun to central body in km;
m_S = 5.683e26;                      % Mass of perturbing body (Saturn) in kg
SOI = d_Scent*(m_central/m_S)^(2/5);          % Central body SOI in km

g0 = 9.81e-3;                        % Gravitational acceleration at Earth's surface in km/s^2

%% User inputs -  Unit vectors
% Note that typically the defaults are fine except the distance unit should
% depend on the scales involved, eg orbital distance from central body

MU = m0;                             % Initial mass in kg
DU = 1e5;                        % Distance unit in km
TU = sqrt((DU^3)/mu_central);        % Time unit in s
VU = DU/TU;                          % Velocity unit in km/s

%% User inputs - Bounds on decision variables

% For initial phase: initial epoch, initial RA DEC, initial vinf
% For each phase: tof, mf, vinfi and vinff for the gravity assist at the
% end of the trajectory, control vector u. Vinff not needed for final
% phase.
% For final phase: vinfi but no vinff because not doing a GA

h_mins = [50];                     % Minimum height above flyby body in km
h_maxs = [1000];                   % Maximum height above flyby body in km

t_ephem_min = 'January 1 2035';      % Time to start fetching planetary ephemeris data, must be string of form Month d yyyy
t_ephem_max = 'June 1 2035';      % End of time ephemeris data taken over

t0_min = 'January 1 2035';           % Earliest trajectory start epoch, string form Month d yyyy
t0_max = 'April 1 2035';           % Latest trajectory start epoch, string form Month d yyyy

% If start on SOI
RAi_min = 0;                         % Min initial RA in rad
RAi_max = 2*pi;                      % Max inital RA in rad
DECi_min = 0;                        % Min initial DEC in rad
DECi_max = 2*pi;                     % Max initial DEC in rad

% If start in orbit
% ai_min = 5000/DU;                    % Min semimajor axis in DU
% ai_max = 1e5/DU;                     % Max semimajor axis in DU
% AoPi_min = -2*pi;                    % Min argument of perigee of initial orbit in rad
% AoPi_max = 2*pi;                     % Max argument of perigee of initial orbit in rad
% ei_min = 0.8;                        % Min eccentricity of initial orbit, unitless
% ei_max = 0.99;                       % Max eccentricity of initial orbit, unitless
% inci_min = 0.680678;                 % Min inclination of initial orbit in rad
% inci_max = 0.680678;                 % Max inclination of initial orbit in rad
% RAANi_min = -2*pi;                   % Min right ascension of ascending node of initial orbit in rad
% RAANi_max = 2*pi;                    % Max right ascension of ascending node of initial orbit in rad
% TAi_min = 0;                         % Min true anomaly in of initial orbit rad
% TAi_max = 2*pi/4;

% If end on SOI
% RAf_min = 0;                          % Min initial RA in rad
% RAf_max = 2*pi;                       % Max inital RA in rad
% DECf_min = 0;                         % Min initial DEC in rad
% DECf_max = 2*pi;                      % Max initial DEC in rad

% Velocity only needs one bound, because the negative of their value is the
% minimum and the positive is the maximum bound. All in VU
vi_bound = 5/VU;                      % Abs value of max/min initial velocity relative to central body in VU
viRel_bound = 2/VU;                   % Initial velocity relative to a flyby body
vf_bound = 2/VU;                      % Abs value of max/min final velocity relative to central body in VU
vfRel_bound = 8/VU;                   % Final velocity relative to a flyby body
vinf_bound = [3,1.5]/VU;                % Velocity bounds for all GA v_inf values

% Bounds on points in free space. These are Cartesian coordinates defined
% relative to the central body. Since space is so big, only use this if you
% know the point in free space you want to reach quite well, eg a point on
% a nominal trajectory.
x_free_min = 1.1e8/DU;
y_free_min = 1.1e8/DU;
z_free_min = 1/DU;
x_free_max = 1.11e8/DU;
y_free_max = 1.11e8/DU;
z_free_max = 1/DU;
vx_free_min = 5/VU;
vy_free_min = 5/VU;
vz_free_min = 5/VU;
vx_free_max = 15/VU;
vy_free_max = 15/VU;
vz_free_max = 15/VU;

% Set bounds for things present in every phase, one value for each per
% phase
dt_min = [1,1]*86400/TU;          % Minimum time of flight in TU
dt_max = [1.5,1.5]*86400/TU;          % Maximum time of flight in TU
mf_min = [12.5,12.5]/MU;                % Minimum final mass in MU
mf_max = [14,14]/MU;                % Maximum final mass in MU
function [r, v] = orb2state(mu, orbels)
% Inputs: mu = gravitational parameter of the central body in km^3/s^2
% orbels = Vector containing the Keplerian orbital elements of the orbit.
% 1. semimajor axis in km, 2. argument of perigee in radians, 3.
% eccentricity, 4. inclination in radians, 5. right ascension of ascending
% node in radians, 6. true anomaly in radians

% Unpack Keplerian elements from input
a = orbels(1);            % Semimajor axis
AoP = orbels(2);           % Argument of periapsis omega in radians
e = orbels(3);             % Eccentricity
inc = orbels(4);             % Inclination in radians
RAAN = orbels(5);          % Right ascension of ascending node Omega in radians
t_anom = orbels(6);         % True anomaly nu in radians

% Find magnitudes
r_mag = (a*(1 - e^2))/(1 + e*cos(t_anom));  % Magnitude of position
v_mag = sqrt(mu/(a*(1-e^2)));

% Position and velocity in the perifocal reference frame
rp = r_mag*[cos(t_anom);
            sin(t_anom);
            0];

vp = v_mag*[-sin(t_anom);
            e + cos(t_anom);
            0];

% Disgusting matrices
% Transformation matrix from argument of perigee
trans_aop = [cos(AoP), sin(AoP), 0;
            -sin(AoP), cos(AoP), 0;
            0,         0,        1];

% Transformation matrix from inclination
trans_i = [1, 0,       0;
           0, cos(inc),  sin(inc);
           0, -sin(inc), cos(inc)];

% Transformation matrix from right ascension of ascending node
trans_RA = [cos(RAAN),  sin(RAAN), 0;
            -sin(RAAN), cos(RAAN), 0;
            0,          0,         1];

% Stick em all together into one gross transformation matrix
trans_mat = (trans_aop*trans_i*trans_RA).';

% Transform the position and velocity from the perifocal frame to the
% geocentric equatorial frame
r = trans_mat*rp;
v = trans_mat*vp;
function orbels = state2orb(mu, rv)
% Inputs: mu = gravitational parameter of the central body in km^3/s^2
% orbels = Vector containing the Keplerian orbital elements of the orbit.
% 1. semimajor axis in km, 2. argument of perigee in radians, 3.
% eccentricity, 4. inclination in radians, 5. right ascension of ascending
% node in radians, 6. true anomaly in radians

x = rv(1);
y = rv(2);
z = rv(3);
vx = rv(4);
vy = rv(5);
vz = rv(6);

r = [x,y,z];
v = [vx,vy,vz];

rmag = norm(r);
vmag = norm(v);

% Angular momentum
h = cross(r,v);
hmag = norm(h);

% Node vector, z component of angular momentum?
nvec = cross([0,0,1],h);
n = norm(nvec);

% Eccentricity
evec = ((vmag^2 - mu/rmag)*r - dot(r,v)*v)/mu;
e = norm(evec);

% Mechanical energy
E = 0.5*vmag^2 - mu/rmag;

% Semimajor axis and p
if e ==1
    a = inf;
    p = (hmag^2)/mu;
else
    a = -mu/(2*E);
    p = a*(1 - e^2);
end

% Other elements
inc = acos(h(3)/hmag);
RAAN = acos(nvec(1)/n);
AoP = acos(dot(nvec,evec)/(n*e));
t_anom = acos(dot(evec,r)/(e*rmag));

% Some conditions
if nvec(2) < 0
    RAAN = 2*pi - RAAN;
end

if evec(3) < 0
    AoP = 2*pi - AoP;
end

if dot(r,v) < 0
    t_anom = 2*pi - t_anom;
end

% Put orbital elements into array
orbels = [a, AoP, e, inc, RAAN, t_anom];
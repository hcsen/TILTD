function rdiff = does_tof_reach_distance(dt, consts)

% Get constants from input
N_perts = consts(1);                % Number of perturbing bodies
N_ephem = consts(2);                % Number of ephemeris data points
N_int = consts(3);                  % Number of timesteps for data in integration
dist = consts(4);                   % Distance aiming for
ri = consts(5:7);                   % Initial position
vi = consts(8:10);                  % Initial velocity
ti = consts(11);                    % Initial epoch in TU
TU = consts(12);                    % Time unit in central body frame
mu_central = consts(13);               % Gravitational parameter of central body in DU^3/TU^2
mu_perts = consts(14:13+N_perts);   % Gravitational parameter of perturbing bodies in DU^3/TU^2

et = consts(14+N_perts : 13+N_perts+N_ephem);   % Ephemeris times in s
stateLong = zeros(3*N_perts, N_ephem);
for i = 1:N_perts
    for j = 1:3
        stateLong(j+(i-1)*3,:) = consts(14+N_perts+(j+3*(i-1))*N_ephem : 13+N_perts+(j+3*(i-1)+1)*N_ephem);
    end
end

%% Spline for ephemeris data over trajectory only

times = linspace(ti-abs(dt), ti+abs(dt), N_int)*TU;
pert_allx = zeros(N_perts, N_int);
pert_ally = zeros(N_perts, N_int);
pert_allz = zeros(N_perts, N_int);
for j = 1:N_perts
    pert_allx(j, :) = spline(et, stateLong(1+(j-1)*3, :), times);
    pert_ally(j, :) = spline(et, stateLong(2+(j-1)*3, :), times);
    pert_allz(j, :) = spline(et, stateLong(3+(j-1)*3, :), times);
end

%% Propagate trajectory by explicit integration for the guessed dt

X0 = [ri, vi];
% dt/86400*TU
tspan = [ti, ti+dt];

[t, state] = ode45(@(t,X0) EoM_hifi_coast(t, X0, mu_central, N_perts, mu_perts, ...
    pert_allx, pert_ally, pert_allz, times, TU), tspan, X0);

% Check if going backwards or forwards in time
if dt < 0
    xf = state(end,1);
    yf = state(end,2);
    zf = state(end,3);
else
    xf = state(1,1);
    yf = state(1,2);
    zf = state(1,3);
end

rf = [xf, yf, zf];

rdiff = norm(rf) - dist;
% figure(2)
% plot3(state(:,1)*1.323808947813119e+04, state(:,2)*1.323808947813119e+04, state(:,3)*1.323808947813119e+04, 'r-')
% 
% figure(4)
% vGA = [state(:,4), state(:,5), state(:,6)];
% plot(t*TU, vecnorm(vGA,2,2)*0.608632678338646, 'r-')
end
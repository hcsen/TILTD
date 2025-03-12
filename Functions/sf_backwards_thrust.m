function [t,r,v,m] = sf_backwards_thrust(rf, vf, dt_p, tf, mf, u, x_pert, y_pert, z_pert, sfconsts)
% % Backwards shooter for the Sims-Flanagan transcription, with two
% perturbing bodies and no thrusting.

%% Get constants and set up matrices from input

% Extract constants from sfconsts
N = sfconsts(1);
mu_central = sfconsts(2);
T = sfconsts(3);
Isp = sfconsts(4);
g0 = sfconsts(5);
N_act = sfconsts(6);
D = sfconsts(7);
N_perts = sfconsts(8);
mu_perts = sfconsts(9:N_perts+8);

% Create empty matrices to store position, velocity, time, mass, and
% thrust. Need a position, velocity, and time for the centre of each
% segment, control node, and match point.
r = zeros(N/2+2, 3);                  % Position relative to the central body in DU
r_mag = zeros(1, N/2+2);         % Distance magnitude to centre of central body in DU
v = zeros(N/2+2, 3);                  % Velocity just before impulse, at node, and match point in VU
v_mag = zeros(1,N/2+2);
t = zeros(1, N/2+2);                   % Epoch in TU
m = zeros(1, N/2+1);
                                                   
% Put initial values in those arrays. Since this is backwards shooting,
% these are the values being targeted at the end of the trajectory.
r(end, :) = rf;                 % Position relative to central body, in DU
v(end, :) = vf;                 % Velocity in VU
t(end) = tf;                     % Final epoch in TU
m(end) = mf;

% Sort out that only using the latter half of perturbing body ephemerides
% and thrust control vector
x_pertUse = zeros(N_perts, N/2);
y_pertUse = zeros(N_perts, N/2);
z_pertUse = zeros(N_perts, N/2);
for i = 1:N_perts
    x_pertUse(i,:) = x_pert(i,N/2+1:N);
    y_pertUse(i,:) = y_pert(i,N/2+1:N);
    z_pertUse(i,:) = z_pert(i,N/2+1:N);
end
u = u(:,N/2+1:N);

%% This loop moves the trajectory backwards through one segment with each
% iteration. Segments are numbered k and go from N/2+1 to 2, because there
% are N segments in the full trajectory and the forwards-shooting part is
% the first half.

for k = (N/2+1):-1:1
    
%% First, define some values for the current position, velocity, and time, for this segment

    % From the array of all positions and velocities, get the position and
    % velocity immediately after the impulse in the previous segment, in DU
    % and VU. For k = N/2+1, these are the initial position and velocity
    r1 = r(k+1, :).';
    v1 = v(k+1, :).';

    % Magnitude of position and velocity immediately after impulse in the
    % previous segment, in DU and VU
    r_mag(k+1) = norm(r1);
    v_mag(k+1) = norm(v1);

    % Length of time for segment to occur, in TU. First and last are
    % half-segments, the rest are just the total phase time divided by the
    % number of segments
    if k == N/2+1 || k == 1
        dt = (dt_p)/(2*N);
    else
        dt = (dt_p)/N;
    end

    % Current epoch in TU - halfway through the segment arc, at the impulse
    t(k) = t(k+1) - dt;
    
    % Alpha = 1/semimajor axis, from the Vis-Viva equation, in DU^-1
    alpha = 2/r_mag(k+1) - (v_mag(k+1)^2)/mu_central;
    
%% Initial guess for universal anomaly chi.
% Chi is used to solve Kepler's equation to propagate position and
% velocity. It is found iteratively, so an initial guess is required. The
% guess made depends on the type of orbit.

    % For the elliptical case, alpha is positive
    if alpha > 1e-12
        % For ellipses, chi = theta/sqrt(alpha), where theta is the
        % angle between the initial and final position vectors.
        chi = alpha * sqrt(mu_central) * dt;

    % Hyperbolic case, where alpha is negative
    elseif alpha < -1e-12
        % In the hyperbolic case, chi = theta/sqrt(-alpha). I don't
        % know where this initial guess comes from!
        chi = (sqrt(mu_central)/(10*r_mag(k+1))) * dt;

    % Parabolic case, where alpha = 0. Found tolerance of +/- 1e-12
    % is needed as it'll never be perfectly 0
    else
        % In the parabolic case, chi = sigma1 - sigma0, where
        % sigma = r.v/sqrt(mu). First guess for chi. Again, not
        % sure where this comes from
        chi = (sqrt(mu_central)/(10*r_mag(k+1))) * dt;
    end                         % End of loop for initial guess for chi

%% Find universal anomaly chi for solving Kepler's equation. Uses
% Laguerre-Conway method developed by Der, shown to converge within 5 iterations

    for i = 1:5
        % Universal variables are defined differently depending on the type
        % of orbit, so check what type of orbit and apply the right
        % equations
        
        % For the elliptical case
        if alpha > 1e-12
            % Functions used to find universal variables. They don't mean
            % anything physically, they just shorten the equations for
            % universal variables.
            y = alpha * chi^2;
            C = (1/y) * (1 - cos(sqrt(y)));
            S = (1/((y^1.5))) * (sqrt(y) - sin(sqrt(y)));

            % Universal variables
            U1 = chi * (1 - y*S);
            U2 = chi^2 * C;
            U3 = chi^3 * S;
            U0 = 1 - alpha*U2;

        % Hyperbolic case
        elseif alpha < -1e-12
            % Universal variables, different to elliptical case
            U0 = cosh(sqrt(-alpha)*chi);
            U1 = 1/(sqrt(-alpha)) * sinh(sqrt(-alpha)*chi);
            U2 = (1/alpha) * (1 - U0);
            U3 = (1/alpha) * (chi - U1);

        % Parabolic case
        else
            % Universal variables, again different
            U0 = 1;
            U1 = chi;
            U2 = 0.5*U1*chi;
            U3 = (1/3)*U2*chi;
        end                     % End of loop to find universal variables
        
%% Use universal variables to update the value of chi this iteration, 
% and to propagate r and Kepler's equation backwards a segment

        % Sigma, used to simplify equations, at initial point (future
        % segment)
        sigma1 = dot(r1, v1)/sqrt(mu_central);
        
        % Propagate sigma to the next step (previous segment)
        sigma0 = sigma1*U0 - (1 - alpha*r_mag(k+1))*U1;
        
        % Propagate magnitude of r to just after previous impulse
        r_mag(k) = r_mag(k+1)*U0 - sigma1*U1 + U2;
        
        % Kepler's equation as a function of chi at the initial point
        f = r_mag(k+1)*U1 - sigma1*U2 + U3 - sqrt(mu_central)*dt;
        
        % First and second derivatives of f with respect to chi
        df = r_mag(k);
        ddf = sigma0;
        
        % eta to simplify finding the amount chi changes each iteration.
        % 5 is there as the number of iterations is 5. Replace 5 if a
        % different number of iterations is used, but should converge in 5.
        % Absolute value is taken to avoid imaginary results, Conway 1986
        % mentioned that isn't a problem
        eta = abs((5 - 1)^2 * df^2 - 5*(5 - 1)*f*ddf);

        % Find change in chi
        dchi = (5*f)/(df + sqrt(eta));

        % Update guess for chi
        chi = chi - dchi;
    end                 % End loop to iteratively find universal anomaly chi
    
%% Propagate position and velocity to just before the impulse using Kepler

    % Propagate magnitude of r to just before the next impulse
    r_mag(k) = r_mag(k+1)*U0 - sigma1*U1 + U2;
    
    % Compute Lagrange F and G coefficients, and their time derivatives
    % (F_dot and G_dot)
    F = 1 - U2/r_mag(k);
    F_dot = -((sqrt(mu_central))/(r_mag(k+1)*r_mag(k))) * U1;
    G = 1/sqrt(mu_central) * (r_mag(k)*U1 + sigma0*U2);
    G_dot = 1 - U2/r_mag(k+1);

    % Propagate r and v forwards one segment to just before the impulse of the current segment
    r0  = G_dot*r1 - G*v1;
    v0 = -F_dot*r1 + F*v1;

    % Magnitude of position (distance to central body) just before impulse
    % at current segment
    r_mag(k) = norm(r0);
    
%% At impulse, apply instantaneous change in velocity due to third-body gravitational perturbations
    
    % No impulse at the match point
    if k ~= 1
        % Distances supplied relative to central body, find distance from
        % spacecraft to perturbing body
        r_pertsc = zeros(3,N_perts);
        r_pert = zeros(3,N_perts);
        pert_SC = zeros(3,N_perts);
        pert_body = zeros(3,N_perts);
        for i = 1:N_perts
            r_pertsc(1,i) = x_pertUse(i,k-1) - r0(1);
            r_pertsc(2,i) = y_pertUse(i,k-1) - r0(2);
            r_pertsc(3,i) = z_pertUse(i,k-1) - r0(3);
            r_pert(1,i) = x_pertUse(i,k-1);
            r_pert(2,i) = y_pertUse(i,k-1);
            r_pert(3,i) = z_pertUse(i,k-1);
            
            % Perturbation on spacecraft
            pert_SC(:,i) = (mu_perts(i)*r_pertsc(:,i))/(norm(r_pertsc(:,i))^3);
            pert_body(:,i) = (mu_perts(i)*r_pert(:,i))/(norm(r_pert(:,i))^3);
        end
        
        v0 = v0 - (sum(pert_SC,2)*dt - sum(pert_body,2)*dt);
        
        % Thrust occurs at impulses too
        dm = T/(Isp*g0);
        dv_max = (N_act * D * T* dt)/(m(k));
        
        % Change in velocity due to thrust
        m(k-1) = m(k) + norm(u(:,k-1)) * D*dt*dm;
        v0 = v0 - dv_max*u(:,k-1);
    end
    
    r(k,:) = r0;
    v(k,:) = v0;
end         % End loop propagating through each segment


end
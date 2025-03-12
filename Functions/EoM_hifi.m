function dXdt = EoM_hifi(t, X0, mu_central, N_perts, mu_perts, x_centPertBig, y_centPertBig, z_centPertBig, times, TU)
% Equations of motion to take a high-fidelity trajectory from an initial
% state vector propagated through a set time of flight, with perturbations
% due to n bodies and thrusting

r = [X0(1), X0(2), X0(3)];
norm_r3 = norm(r)^3;

dv_pert = zeros(3, N_perts);

% Perturbing velocity due to each perturbing body
for i = 1:N_perts
    % Interpolate to get perturbing body relative to central body at the current time
    x_CP = interp1(times, x_centPertBig(i,:), t*TU);
    y_CP = interp1(times, y_centPertBig(i,:), t*TU);
    z_CP = interp1(times, z_centPertBig(i,:), t*TU);
    r_CP = [x_CP, y_CP, z_CP];
    
    % Distance from spacecraft to perturbing body
    r_pert = r_CP.' + r;
    rp_norm3 = norm(r_pert)^3;
    
    % Work out the change in velocity due to perturbing body
    dv_pert(:,i) = (mu_perts(i)*r_pert)/rp_norm3 - (mu_perts(i)*r_CP.')/(norm(r_CP)^3);
end

dXdt(1) = X0(4);
dXdt(2) = X0(5);
dXdt(3) = X0(6);
dXdt(4) = -(mu_central*r(1))/norm_r3 + sum(dv_pert(1,:));
dXdt(5) = -(mu_central*r(2))/norm_r3 + sum(dv_pert(2,:));
dXdt(6) = -(mu_central*r(3))/norm_r3 + sum(dv_pert(3,:));
dXdt = dXdt.';
end
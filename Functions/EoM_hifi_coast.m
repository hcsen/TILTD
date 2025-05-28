function dXdt = EoM_hifi_coast(t, X0, mu_central, N_perts, mu_perts, x_centPertBig, y_centPertBig, z_centPertBig, times, TU)
% Equations of motion to take a high-fidelity trajectory from an initial
% state vector propagated through a set time of flight, with perturbations
% due to n bodies and thrusting

r = [X0(1), X0(2), X0(3)];
norm_r3 = norm(r)^3;

dv_pert = zeros(3, N_perts);
ty = t*TU;
lower = find(times <= ty, 1, 'last');
higher = find(times >= ty, 1);
if ty > times(end)
    higher = length(times);
    lower = length(times)-1;
elseif ty < times(1)
    lower = 1;
    higher = 2;
elseif ty == times(higher) || ty == times(lower)
    if lower ~= 1
        lower = lower - 1;
    elseif higher ~= length(times)
        higher = higher + 1;
    end
end
t1 = times(higher);
t0 = times(lower);

% Perturbing velocity due to each perturbing body
for i = 1:N_perts
    % Interpolate to get perturbing body relative to central body at the current time
%     x_CP = interp1(times, x_centPertBig(i,:), t*TU);
%     y_CP = interp1(times, y_centPertBig(i,:), t*TU);
%     z_CP = interp1(times, z_centPertBig(i,:), t*TU);
    xl = x_centPertBig(i,lower);
    xu = x_centPertBig(i,higher);
    x_CP = xl + (t1 - t0)*((xu - xl)/(t1 - t0));
    yl = y_centPertBig(i,lower);
    yu = y_centPertBig(i,higher);
    y_CP = yl + (t1 - t0)*((yu - yl)/(t1 - t0));
    zl = z_centPertBig(i,lower);
    zu = z_centPertBig(i,higher);
    z_CP = zl + (t1 - t0)*((zu - zl)/(t1 - t0));
    r_CP = [x_CP, y_CP, z_CP];
    
    % Distance from spacecraft to perturbing body
    r_pert = r_CP + r;
    rp_norm3 = norm(r_pert)^3;
    % Work out the change in velocity due to perturbing body
    dv_pert(:,i) = (mu_perts(i)*r_pert)/rp_norm3 - (mu_perts(i)*r_CP)/(norm(r_CP)^3);
end

dXdt(1) = X0(4);
dXdt(2) = X0(5);
dXdt(3) = X0(6);
dXdt(4) = -(mu_central*r(1))/norm_r3 + sum(dv_pert(1,:));
dXdt(5) = -(mu_central*r(2))/norm_r3 + sum(dv_pert(2,:));
dXdt(6) = -(mu_central*r(3))/norm_r3 + sum(dv_pert(3,:));
dXdt = dXdt.';
end
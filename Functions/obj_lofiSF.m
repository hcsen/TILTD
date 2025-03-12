function [f] = obj_lofiSF(x,N)
%% Objective function for MGALT models where
% Minimise tof because there's no mass

mf = x(N);
f = -mf;
end
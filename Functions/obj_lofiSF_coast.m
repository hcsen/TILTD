function [f] = obj_lofiSF_coast(x,N)
%% Objective function for MGALT models where
% Minimise tof because there's no mass

f = x(N);
end
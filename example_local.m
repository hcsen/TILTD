% Example local interactive usage

% Add micepath. (Ideally this should be done outside the function, that way we can remove micepath as a config value)
addpath('C:\Users\hcsen\Documents\MATLAB\mice\src\mice', 'C:\Users\hcsen\Documents\MATLAB\mice');
output = lofi_search('Inputs/verification_conf_base', 'local_conf');

plot_lofiSF(output.best,output.consts);
%% Save result

% Saves the optimised result, constants, constraint violation, lower and
% upper bounds. Can replace with anything else you want to save
% dlmwrite('FILE_NAME_HERE.csv', optimised, 'delimiter', ',', 'precision', 20);
% dlmwrite('FILE_NAME_HERE_consts.csv', consts, 'delimiter', ',', 'precision', 20);
%dlmwrite('viol.csv', violation_archive, 'delimiter', ',', 'precision', 20); % This is temporary, for purpose of verification.
% dlmwrite('FILE_NAME_HERE_lb.csv', lb, 'delimiter', ',', 'precision', 20);
% dlmwrite('FILE_NAME_HERE_ub.csv', ub, 'delimiter', ',', 'precision', 20);

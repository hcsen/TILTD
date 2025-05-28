% Example non interactive usage


% Add micepath. (Ideally this should be done outside the function, that way we can remove micepath as a config value)
addpath('../mice/src/mice', '../mice/lib');
output = lofi_search('Inputs/verification_conf_base', 'hpc_conf');
save('output.mat', "output", '-mat');


% Then, interactively call.
% load('output.mat', "output");
% plot_lofiSF(output.best,output.consts);
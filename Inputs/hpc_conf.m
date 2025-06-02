rng_seed = str2num(getenv('SLURM_ARRAY_TASK_ID')); % Needs to be 'shuffle', or int
mice_path = {'..', 'mice'};   % Path to mice library.
kernel_path = {'Kernels', 'saturn_ev.tm'};  % Path to kernel

use_parallel = true;   % whether to use parpool or not.
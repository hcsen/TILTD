#!/bin/bash -e

#SBATCH --account    nesi99991      # Your project code
#SBATCH --job-name   TILDT          # Name to appear in squeue 
#SBATCH --time       01:00:00       # Max walltime 
#SBATCH --cpus-per-task 4
#SBATCH --mem-per-cpu   2G             # Max memory
#SBATCH --partition  milan

module load MATLAB/2023b

# Run the MATLAB script MATLAB_job.m 
matlab -batch "addpath('Inputs'); lofi_search TEE_basic hpc_conf"

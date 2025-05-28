#!/bin/bash -e

#SBATCH --account    nesi99999      # Your project code
#SBATCH --job-name   TILDT          # Name to appear in squeue 
#SBATCH --time       01:00:00       # Max walltime 
#SBATCH --cpus-per-task 40
#SBATCH --mem-per-cpu   3G             # Max memory
#SBATCH --partition  milan

module load MATLAB/2023b

# Run the MATLAB script MATLAB_job.m 
matlab -batch "example_sbatch"

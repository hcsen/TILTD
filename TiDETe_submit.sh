#!/bin/bash -e

##SBATCH --account    uoa04348                    # Your project code
#SBATCH --job-name   TILTD_TiDETe                # Name to appear in queue
#SBATCH --output     %x_%A_%a.out                # Output filename
#SBATCH --time       12:00:00                    # Max walltime 
#SBATCH --cpus-per-task 20
#SBATCH --mem-per-cpu   3G                       # Max memory
#SBATCH --partition  milan
#SBATCH --profile    all
##SBATCH --mail-user  hsen495@aucklanduni.ac.nz   #Email user when job is completed or fails
##SBATCH --mail-type  END,FAIL
#SBATCH --array      1-5                         #Submit job 5 times. Note, doesn't have to start at one.

module purge

module load MATLAB/2023b

export LANG=en_NZ.UTF-8

# Run the MATLAB script MATLAB_job.m 
matlab -batch "TiDETe_run1"

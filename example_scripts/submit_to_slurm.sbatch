#!/bin/bash

#SBATCH --partition=cpu-galvani         # Partition (job queue)
#SBATCH --time=08:00:00            		# Runtime in D-HH:MM
#SBATCH --ntasks=1                		# Number of tasks (set for MPI, for OpenMP to 1)
#SBATCH --job-name=clida	            # Job name
#SBATCH --output=%j.out  				# File to which STDOUT will be written
#SBATCH --error=%j.err   				# File to which STDERR will be written

source $HOME/.bashrc
conda activate clidaEnv

python download_era5.py



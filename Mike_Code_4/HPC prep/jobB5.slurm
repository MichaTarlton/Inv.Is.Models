#!/bin/bash
#SBATCH -J foris5b5  # Sensible name for the job
#SBATCH -N 2             # Allocate 2 nodes for the job
#SBATCH --ntasks-per-node=16
#SBATCH -t 3-00:00:00    # Upper time limit for the job (DD-HH:MM:SS)
#SBATCH -p CPUQ
#SBATCH --mem 25000

module load MATLAB/2019b

matlab -nodisplay -nodesktop -nosplash -r "foris5b5"

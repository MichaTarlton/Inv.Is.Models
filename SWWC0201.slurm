#!/bin/bash
#SBATCH -J SWWC0201  # Sensible name for the job
#SBATCH -N 1             # Allocate nodes for the job
##SBATCH --ntasks-per-node=18 # MPI ONLY!
#SBATCH --account=share/mh-kin
#SBATCH -t 02-00:00:00    # Upper time limit for the job (DD-HH:MM:SS)
#SBATCH -p CPUQ			# setting to GPU because it still has a cpu
#SBATCH --mem=27G # Sey memory allocation. Set to 110g to secure a dedicated node.
#SBATCH -c 20 ## Number of cores

##SBATCH --output=./output/$SLURM_JOB_ID/output.slurm-%A_%a.out

##module load GCC/9.3.0  OpenMPI/4.0.3 HDF5/1.10.6

SBATCH --array=1-132
echo $SLURM_ARRAY_JOB_ID 
echo $SLURM_ARRAY_TASK_ID 
echo $PWD
SCRATCH_DIRECTORY=/home/michaeta/$SLURM_ARRAY_JOB_ID 
mkdir -p $SCRATCH_DIRECTORY
echo $SCRATCH_DIRECTORY


module load MATLAB/2020b 


matlab -nodisplay -nodesktop -nosplash -r "fi0201($SLURM_ARRAY_JOB_ID,$SLURM_ARRAY_TASK_ID)"   ## Can't use -r if passing input into matlab

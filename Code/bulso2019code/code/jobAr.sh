#!/bin/bash -l
#PBS -l walltime=40:00:00
#PBS -l nodes=1:ppn=1

#! General jobscript. Do not modify below
#!echo $SLURM_ARRAY_TASK_ID
#!CHANGE="./${TEXT}${SLURM_ARRAY_TASK_ID}"
#!cd $CHANGE
#!echo $CHANGE
#! NODES='wc -l <$PBS_NODEFILE'
module load matlab
matlab -nodisplay -r "$MLS(${SLURM_ARRAY_TASK_ID});"
module unload matlab



#!cp $DATASET $SCRATCH
#!cd $SCRATCH
#!chkfile $OUTFILE
#!YourProgram $DATASET > $OUTFILE


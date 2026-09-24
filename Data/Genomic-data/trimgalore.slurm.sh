#!/bin/bash
#
#SBATCH -J tg-SAMPLE_ID
#SBATCH -N 1
#SBATCH --partition=hpc
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --ntasks-per-core=2
#SBATCH --mem=10G
#SBATCH --time 10:00:00

#SBATCH --output slurm-%x-%j.out
#SBATCH --error slurm-%x-%j.err

############
# specify variables if needed
SAMPLE_ID=
forward=../rawdata/${SAMPLE_ID}/${SAMPLE_ID}-1.fq.gz
reverse=../rawdata/${SAMPLE_ID}/${SAMPLE_ID}-2.fq.gz

#From here on out the magic happens
echo "Host: $HOSTNAME"
echo "Job: $JOB_ID"

#get going
echo -e "\n[$(date)] - Starting ..\n"

singularity exec -B /cl_tmp docker://chrishah/trim_galore:0.6.0 \
trim_galore \
--paired --length 70 -r1 71 -r2 71 --retain_unpaired --stringency 2 --quality 30 \
$forward $reverse

#notify when done
echo -e "\n[$(date)] - Done!\n"


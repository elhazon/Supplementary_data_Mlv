#!/bin/bash
#
#SBATCH -J dup-SAMPLE_ID
#SBATCH -N 1
#SBATCH --partition=hpc
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --ntasks-per-core=2
#SBATCH --mem=40G
#SBATCH --time 10:00:00

#SBATCH --output slurm-%x-%j.out
#SBATCH --error slurm-%x-%j.err

############
# specify variables if needed
SAMPLE_ID=
bam=#/${SAMPLE_ID}.sorted.Aligned.sortedByCoord.out.bam


#From here on out the magic happens
echo "Host: $HOSTNAME"
echo "Job: $JOB_ID"

#get going
echo -e "\n[$(date)] - Starting ..\n"

##  MARK DUPLICATES  ##
singularity exec -B #/rnaseq-mini_latest.sif \
	picard MarkDuplicates \
	I=${bam} \
	O=${SAMPLE_ID}.markdup.bam \
	M=${SAMPLE_ID}.metrics.markdup.txt \
	REMOVE_DUPLICATES=false VALIDATION_STRINGENCY=SILENT ASSUME_SORTED=true

## REMOVE DUPLICATES ##
singularity exec -B #/rnaseq-mini_latest.sif \
	picard MarkDuplicates \
	I=${bam} \
	O=${SAMPLE_ID}.rmdup.bam \
	M=${SAMPLE_ID}.metrics.rmdup.txt \
	REMOVE_DUPLICATES=true VALIDATION_STRINGENCY=SILENT ASSUME_SORTED=true

#notify when done
echo -e "\n[$(date)] - Done!\n"

#!/bin/bash
#
#SBATCH -J fc-SAMPLE_ID
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
bam_prefix=#/${SAMPLE_ID}
gtf=#/lplatythorax_final_scaffolds_v1_braker3b.gff3

#From here on out the magic happens
echo "Host: $HOSTNAME"
echo "Job: $JOB_ID"

#get going
echo -e "\n[$(date)] - Starting ..\n"


singularity exec -B #PATH/rnaseq-mini_latest.sif \
	featureCounts \
	-T ${SLURM_CPUS_PER_TASK} -p -F GTF -t exon -g Parent \
	-a $gtf \
	-o ${SAMPLE_ID}.markdup.featurecount \
	${bam_prefix}.markdup.bam

singularity exec -B #PATH/rnaseq-mini_latest.sif \
	featureCounts \
	-T ${SLURM_CPUS_PER_TASK} -p -F GTF -t exon -g Parent \
	-a $gtf \
	-o ${SAMPLE_ID}.rmdup.featurecount \
	${bam_prefix}.rmdup.bam

#notify when done
echo -e "\n[$(date)] - Done!\n"

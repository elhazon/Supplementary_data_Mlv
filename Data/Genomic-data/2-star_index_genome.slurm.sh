#!/bin/bash
#
#SBATCH -J STARindex
#SBATCH -N 1
#SBATCH --partition=hpc
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --ntasks-per-core=2
#SBATCH --mem=100G
#SBATCH --time 48:00:00

#SBATCH --output slurm-%x-%j.out
#SBATCH --error slurm-%x-%j.err

############
# specify variables if needed
fasta=#/lplatythorax.final_scaffolds_v1.fasta
gtf=/#/lplatythorax_final_scaffolds_v1_braker3b.gff3
refdir=#
overhang=149 # this should be read length -1 

#From here on out the magic happens
echo "Host: $HOSTNAME"
echo "Job: $JOB_ID"

#get going
echo -e "\n[$(date)] - Starting ..\n"

## note in this case I am calling a local copy of the image 'docker://passdan/rnaseq-mini', because I haven't made it and there is no guarantee that it will still exist
## in the same form in the cloud tomorrow 

singularity exec -B /#PATH/rnaseq-mini_latest.sif \
	STAR \
	--runThreadN ${SLURM_CPUS_PER_TASK} \
        --limitGenomeGenerateRAM $(( SLURM_MEM_PER_NODE * 1000000 )) \
	--runMode genomeGenerate \
	--genomeDir ${refdir} \
	--genomeFastaFiles ${fasta} \
	--sjdbGTFfile ${gtf} \
	--sjdbGTFtagExonParentTranscript Parent \ 
	--sjdbOverhang ${overhang}

#notify when done
echo -e "\n[$(date)] - Done!\n"

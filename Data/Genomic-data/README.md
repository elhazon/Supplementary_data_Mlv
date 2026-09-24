The data here contains the original representation of the genome of _Lasius platythorax_ (now deposited under NCBI GCA_964030505.1 with new sequence headers), including the results of strucutural annotation and gene coding sequences (presented in Feldmeyer et al. 2024. doi: https://doi.org/10.1101/2023.07.18.549505), which was kindly provided by the first author Barbara Feldmeyer. 

We deposit it here to provide full reproducibility of our analyses.

Genome assembly is split up into multiple parts because of file size limitations on Github. To reconstruct the original file:
```bash
zcat lplatythorax.final_scaffolds_v1.00* > lplatythorax.final_scaffolds_v1.fasta
```

Additional files:
- lplatythorax_final_scaffolds_v1_braker3b.gff3.gz - gene coordinates (as predicted by Feldmeyer et al. 2024)
- lplatythorax_final_scaffolds_v1_braker3b.codingseq.gz - coding sequences (as predicted by Feldmeyer et al. 2024)
- lplatythorax_final_scaffolds_v1_braker3b.aa.gz - amino acid sequences (as predicted by Feldmeyer et al. 2024)
- eggnog_results.emapper.annotations.gz - functional annotation of genes obtained with eggNOG 1.0.3 (our study)
- eggnog.GOterms.tsv.gz - GOterms per gene (our study)
- trimgalore.slurm.sh - script for trimming the raw sequencing reads
- 2-star_index_genome.slurm.sh - script for mapping against reference genome
- 4-markduplicates.slurm.sh - script for removing duplicates
- 5-featurecountr.slurm.sh - script for read counts per read



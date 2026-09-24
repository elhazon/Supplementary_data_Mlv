# r-script_Malvolio
Script and raw data files for "Infection induced changes in gene expression linked to foraging choices in self-medicating Lasius platythorax ants"

This repository includes all necessary files and R-script for analysis.
The R-script can be found in the "Script" folder ("full-analysis-script.R").  
The "Data" folder includes all the data files for the analysis.  
The "Genomic-data" folder includes the genome files used for the analysis, as well as scripts.
For the script to run as is, add the script to the same folder as the data files.

# file descriptions

## Script
"full-analysis-script.R"
Contains a commented R-script to fully recreate the original analysis with explanation to each step in detail.

## Data

### Genome-data
Includes all the genomic data used in the analysis, along with scripts.  
The folder includes a separate README.md file with detailed information on the contents.  
___  

### "Dataset1.xlsx" - Contains the raw foraging and survival data for experiment a. 

Sheet "Foraging" column variables:

OC (Original Colony) = The identifier of the original collected wild nest, further split into the experimental colonies.

Box = Experimental colony ID.

Day = Day of experiment.

Time = Time of day (10, 14, 18, 22).  
___  

Sheet "Survival" column variables:

Colony = Experimental colony ID.

Survival = Survival event (1 = dead, 0 = alive).

Day = Day of experiment.

OC = original collected wild nest ID.

Treatment = Infection treatment (C = sham-treated control, F = exposed to pathogen).

Diet = Diet treatments (1 = Standard food, 2 = Food choice, 3 = ROS food).  
___  

### "Dataset2.xlsx" Contains the RNA sequencing data for the differential gene expression and GO-Term analysis (experiment a).

Sheet "lPlatAAseq" contains the full amino acid sequences for each of the genes from the Lasius platythorax genome.

name = Gene name.

sequence = Full amino acid sequence.  
___  

Sheet "sample_info_d2" contains the formatted sample information for the DESeq2 analysis.

SampleName = Sample name.

Treatment = Infection treatment.

Time = Sampling event (1 = 2 days post exposure).

Replicate = Blocking factor for original wild nest ID.  
___  

Sheet "sample_info_d4" contains the formatted sample information for the DESeq2 analysis.

SampleName = Sample name.

Treatment = Infection treatment.

Time = Sampling event (2 = 4 days post exposure).

Replicate = Blocking factor for original wild nest ID.  
___  

Sheet "reads_d2" contains the raw reads from the RNA sequencing for the day 2 samples.

First column (nameless) = Gene name.

Columns JR11 - JR30 = Reads for each gene of each sample.  
___  

Sheet "reads_d4" contains the raw reads from the RNA sequencing for the day 4 samples.

First column (nameless) = Gene name.

Columns JR31 - JR50 = Reads for each gene of each sample.  
___  

Sheet "Freq_tbl_ancestor" contains the data on which GO term is associated with parent processes for the biological processes enriched with upregulated genes on day 4.

GO.Term = GO term.

Parent.Process = Parent process higher in the ontology hierarchy.

fr = a supporting column for addition and process for visualization.  
___  

### "Dataset3.xlsx" contains the raw foraging and gene expression data for the analysis of experiment b.

Sheet "Mlv_foraging" - Contains the raw foraging  data for experiment b. 
Column variables:

OC (Original Colony) = The identifier of the original collected wild nest, further split into the experimental colonies.

Box = Experimental colony ID.

Day = Day of experiment.

Time = Time of day (9, 12, 15, 18, 21).

Foragers = Number of ants in contact with the food.

Food_Type = Types of food in the experiment (1 = standard food, 2 = ROS food).

Treatment = Infection treatment (C = sham-treated control, F = exposed to pathogen).  
___  

Sheet "MCMC" - Contains the raw cq reads for the RT-qPCR analysis.
Column variables:

Colony = Colony ID.

Nest = Original colony ID.

Day = Sampling day.

Treatment = Infection treatment (C = sham-treated control, F = exposed to pathogen).

Sample = Unique sample identifier, shared between replicates.

Plate = PCR plate ID.

MalG = Cq values for the MalG primer reactions.

Act = Cq values for the Actin primer reactions.  
___  

Sheet "meff" - contains the primer efficiencies of the MalG and Act primers.
Column variables: 

gene = gene name.

efficiency = primer efficiency.  
___  

### "combined.byGO.nr.tsv" Contains each gene in the analysis and its associated GO-terms.  
___

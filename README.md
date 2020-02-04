# Supplementary Note from NRG review  

As part of this review, we performed limited meta-analysis with long-read datasets as part of the comparisons of read length, accuracy, and homopolymer analysis. All datasets are publicly available, and we describe briefly how the data were analyzed using the code within this repository. 

## Data Sources
Pacific Biosciences (PacBio) HG002 CLR data was retrieved from PacBio’s public database (https://github.com/PacificBiosciences/DevNet/wiki/HG002-Structural-Variant-Analysis-with-CLR-data) at the following URL: https://downloads.pacbcloud.com/public/dataset/SV-HG002-CLR/m64013_190124_221354.subreads.bam. Library preparation was performed with the SMRTbell Express 2.0 kit and size-selected to be greater than 30 kbp with a BluePippin instrument. 

PacBio CHM13 HiFi data1 was retrieved from SRA accessions SRR9087597, SRR9087598, SRR9087599, and SRR9087600.

Oxford Nanopore Technologies (ONT) CHM13 whole-genome sequencing data was retrieved from the T2T consortium’s GitHub repository (https://github.com/nanopore-wgs-consortium/chm13) at the following URL: https://s3.amazonaws.com/nanopore-human-wgs/chm13/nanopore/rel3/rel3.fastq.gz. This dataset is comprised of long reads (N50 = 35.2 kbp) generated on the PromethION at the University of California, Davis and ultra-long reads (N50 = 146.1 kbp) generated on the MinION and GridION at the University of Washington. The entire dataset was base called with Guppy 3.1.5 using the flipflop  model and separated into long- and ultra-long -read datasets based on the read IDs provided from each institution, which are available at the following URLs: https://s3.amazonaws.com/nanopore-human-wgs/chm13/nanopore/rel3/ids/uwashington.ids.gz and https://s3.amazonaws.com/nanopore-human-wgs/chm13/nanopore/rel3/ids/ucd.ids.gz. 

## Estimation of read length, accuracy, and homopolymer error
PacBio CLR and HiFi data were  aligned to GRCh38 using pbmm2 (https://github.com/PacificBiosciences/pbmm2) v1.1.0 with the following parameters: `--preset SUBREAD`. ONT long- and ultra-long -read data was were aligned to GRCh38 using minimap22 v2.17 with the following parameters: `-ax map-ont -L --eqx`. All data was filtered to exclude secondary, supplementary, and unmapped alignments using SAMtools3 v1.9 and SAM flag 2308. 

Read length was determined from the sequence length column in the SAM file for aligned reads, and read accuracy was calculated from the CIGAR string in the SAM file using the following formula: 100 * (number of matching bases)/(number of matching bases + number of mismatched bases + number of bases in insertions + number of bases in deletions). The code to reproduce these results is available in the GitHub repository listed below in the snakemake `length_and_qv.smk`.

Differences in homopolymer lengths were calculated by first tabulating the homopolymers in the reference genome assembly using the script `Find_homopolymers.py` and then comparing the length of these homopolymers to the aligned reads using the snakemake `homopolymer.smk`. Both programs are available in the GitHub repository listed below. Specifically, the CIGAR operations over homopolymers were used to detect inaccuracies in length and then encoded into a matrix of counts where the row value was the reference genome assembly homopolymer length and the column value was the observed read homopolymer length. 

The datasets generated in the analyses above were then visualized using Matplotlib4 v2.1.2 and seaborn5 v0.10.0 in a Jupyter notebook called `plots_notebook.ipynb` to create Figure 3c,d. 

These analyses were also used to calculate the read lengths, accuracies, and homopolymer profiles in Figure S1 with two exceptions: 1) the reads were aligned to the CHM13 genome assembly v0.6 from Ref. 6, available at the following URL: https://s3.amazonaws.com/nanopore-human-wgs/chm13/assemblies/chm13.draft_v0.6.fasta.gz, and 2) accuracy and homopolymer estimations were restricted to the X chromosome due to the high level of curation performed on this chromosome. 



snakemake -j 160 -p -s length_and_qv.smk --restart-times 1

./Find_Homopolymers.py ~/assemblies/hg38/ucsc.hg38.no_alts.fasta  --chr chrX --out results/homopolymer.smk
snakemake -j 160 -p -s homopolymer.smk --restart-times 3

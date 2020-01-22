import pandas as pd
import numpy as np
import pysam


#
# INPUTS, change as needed
#

# read in / set up inputs
MAX_HP=100 # maximum detected homopolymer
configfile: "config.yaml"
bams = config["bams"]
fai = config["fai"]

print("bams to use")
print(bams)


#
# wildcards
#
contigs = pd.read_csv(fai, names=["contig", "length", "x", "y", "z"] , sep="\t")
contigs = list(contigs["contig"])
IDs = list(range(len(contigs)))
techs = list(bams.keys())

wildcard_constraints:
	ID = "\d+",
	TECH = "|".join(techs)


#
# input functions
#
def get_bam(wildcards):
	return(bams[ str(wildcards.TECH) ])

def get_contig(wildcards):
	return(contigs[int(str(wildcards.ID))])

def get_lengths(bam, contig):
	EQ=7;X=8;I=1;D=2
	bam = pysam.AlignmentFile(bam)
	lengths = []
	for idx, read in enumerate(bam.fetch(contig=contig)):
		counts, events = read.get_cigar_stats()
		if(counts[EQ] > 0 ):
			matchID = counts[EQ]/(counts[EQ] +counts[X])
			indelID = counts[EQ]/(counts[EQ] +counts[X] + counts[I] + counts[D])
			lengths.append((len(read.seq), read.query_name, matchID, indelID))
		else:
			print(read.query_name, len(read.seq), counts, read.cigarstring)
			exit()
		if(idx %100 == 0): 
			sys.stderr.write("\r{}".format(idx))
	bam.close()
	return(lengths)

#
# rules 
#
rule all:
	input:
		txt = "results/all.txt",
		txts = expand("results/{TECH}.txt",  TECH=techs),

rule lengths:
	input:
		bam = get_bam,
	output:
		txt = temp("temp/{TECH}_{ID}.txt"),
	params:
		rgn = get_contig,
	resources:
		mem=lambda wildcards, attempt: 4 + 4*attempt,
	run:
		lengths = get_lengths( input["bam"], params["rgn"] )
		f = open(output["txt"], "w+")
		for read_l, read, matchID, indelID in lengths:
			f.write("{}\t{}\t{}\t{}\t{}\t{}\n".format(read_l, read, matchID, indelID, wildcards.TECH , params["rgn"]))
		f.close()

rule merge1:
	input:
		txt = expand("temp/{{TECH}}_{ID}.txt",  ID=IDs),
	output:
		txt = "results/{TECH}.txt"
	resources:
		mem=32, 
	shell:"""
cat {input.txt} > {output.txt}
"""

rule merge2:
	input:
		txt = expand("results/{TECH}.txt",  TECH=techs),
	output:
		txt = "results/all.txt",
	resources:
		mem=32,
	shell:"""
cat {input.txt} > {output.txt} 
"""



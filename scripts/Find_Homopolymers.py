#!/usr/bin/env python
import argparse
import sys

parser = argparse.ArgumentParser(description="")
parser.add_argument("infile", help="positional input")
parser.add_argument("-n", "--number", help="min homopolymer length", type=int, default=2)
parser.add_argument("-o", "--out", help="pkl with homopoyler regions", default="homopolymers.bed.pkl" )
parser.add_argument("-c", "--chr", help="chromosome to calculate homopolymers on", default="chrX" )
parser.add_argument('-d', help="store args.d as true if -d",  action="store_true", default=False)
args = parser.parse_args()

import pandas as pd
from Bio import SeqIO

min_homo = args.number
homos = {"contig":[], "start":[], "end":[]}

ref = SeqIO.parse(args.infile, "fasta")

for rec in ref:
	if(rec.id != args.chr): continue 
	pre = ""
	start = 0
	idx = start
	length = 1
	seq = str(rec.seq).upper()	
	total = len(seq)
	while idx < total:
		cur = seq[idx]
		if(cur == pre and cur != "N"):
			length += 1
		else:
			if(length >= min_homo):
				homos["contig"].append(rec.id)
				homos["start"].append(idx-length)
				homos["end"].append(idx)
			length = 1
	

		if(idx % 10000 == 0):
			sys.stderr.write("\r{}\t{:.2f}%".format(rec.id, idx/total*100))

		idx += 1
		pre = cur
	sys.stderr.write("\n")

tmp=pd.DataFrame(homos)
tmp.to_pickle(args.out)


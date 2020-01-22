import pysam
import sys
import os
import pandas as pd
import numpy as np

# read in / set up inputs
MAX_HP=100 # maximum detected homopolymer 
configfile: "config.yaml"
bams = config["bams"]
f_homos = config["hp_pkl"]

print("bams to use")
print(bams)


df = pd.read_pickle(f_homos)
rgns = list(zip(df.contig, df.start, df.end))
sys.stderr.write("read in pkl\n")

def batch_it(a, n):
	k, m = divmod(len(a), n)
	return list( (a[i * k + min(i, m):(i + 1) * k + min(i + 1, m)] for i in range(n)) )

n = min(500, len(rgns)) 
batches = batch_it(rgns, n)
IDs = list(range(n))

sys.stderr.write("batches made\n")


techs = list(bams.keys())



def get_bam(wildcards):
	return(bams[ str(wildcards.TECH) ])

wildcard_constraints:
	ID = "\d+",
	TECH = "|".join(techs)


# updated calling
def read_bam_for_homo(bam_f, homos, saveto):
    bam = pysam.AlignmentFile(bam_f)
    chrm, start, end = (homos[0][0], homos[0][1], homos[-1][2])
    sys.stderr.write("{}\t{}\t{}\n".format(chrm, start, end))
    
    data  = { "chrm":[], "start":[], "end":[], "observations":[] }
    
    homo_pos = 0
    cur_homo = homos[homo_pos]
    reads = {}
    total = end-start 
    cur_count = 0

    for pileupcolumn in bam.pileup(chrm, start, end,
		 truncate=True, stepper="nofilter", max_depth=50000, 
			ignore_orphans=False, min_base_quality=0, min_mapping_quality=0):
        # get the position we are at in the chrm
        pos = pileupcolumn.pos 
        #if(pos % 100 == 0): sys.stderr.write("\r{}\t{}".format(chrm, pos))
        # progress through the list of homopolmers until its end is greater than our current position
        while( pos >= cur_homo[2] ): 
            homo_pos += 1
            cur_homo= homos[homo_pos]
            reads = {}
            
            
        if( pos >= cur_homo[1] and pos < cur_homo[2] ): # we are in a homopoymer 
            #sys.stderr.write("{}\t{}\t{}\n".format(pos, cur_homo, reads))
            
            # count where the read has an indel and add it to the read dict
            for pileupread in pileupcolumn.pileups:
                rname = pileupread.alignment.query_name
                if(rname not in reads): reads[rname] = {"count":0, "diff":0 }
                reads[rname]["count"] += 1
                indel = pileupread.indel
                reads[rname]["diff"] += indel
                #sys.stderr.write("{}\n".format(reads))
                
            # if we are the the end of a homopolmer save the results
            if(pos == cur_homo[2] - 1 ):
                length = cur_homo[2] - cur_homo[1]
                lengths = [] 
                for key in reads:
                    #sys.stderr.write("\r{}".format(key))
                    if( length == reads[key]["count"]):
                        diff = reads[key]["diff"]
                        #sys.stderr.write("\r{}".format(diff))
                        lengths.append(int(length + diff))
                data["chrm"].append(cur_homo[0])
                data["start"].append(cur_homo[1])
                data["end"].append(cur_homo[2])
                data["observations"].append(lengths)
                #sys.stderr.write("\r{}\t{}\t{}\t{}".format( cur_homo[0], cur_homo[1], cur_homo[2], len(lengths) ) )

        cur_count +=1 
        if(cur_count % 1000 ==0 ):			
            sys.stderr.write("\r{:.02%}".format( cur_count/total ) )
            
    pd.DataFrame(data).to_pickle(saveto)
    sys.stderr.write("\nWritten {}\n".format(saveto))
    bam.close()

# makes a matrix, where we have counts for observed homoopolmers vs theri real length  
def make_mat(df, max_hp=MAX_HP, extend=10):
    data = np.zeros((max_hp, max_hp+extend))
    df["length"] = df.end-df.start
    for idx, row in df.iterrows():
        sys.stderr.write("\r{}".format(idx))
        for obs in row["observations"]:
            i = row["length"]
            if( i >= 0 and i < max_hp and obs >=0 and obs < max_hp+extend):
                data[row["length"]][obs] += 1
    return(data)


rule all:
	input:
		tbls = expand("results/{TECH}.pkl", TECH=techs),
		mats = expand("results/{TECH}.mat.npy", TECH=techs),

rule get_data:
	input:
		bam = get_bam,
		in_file = f_homos, 
	output:
		pkl = temp("temp/{TECH}_{ID}.pkl"),
	resources:
		mem=lambda wildcards, attempt: 4 + 4*attempt,
	run:
		ID = int(str(wildcards.ID))
		homos = batches[ID]
		read_bam_for_homo(input["bam"], homos, output["pkl"])


rule make_mat:
	input:
		pkl = "temp/{TECH}_{ID}.pkl",
	output:
		npy = temp("temp/{TECH}_{ID}.mat.npy"),
	resources:
		mem=lambda wildcards, attempt: 4 + 4*attempt,
	run:
		df = pd.read_pickle(input["pkl"])
		mat  = make_mat(df)
		np.save(output["npy"], mat)

rule merge:
	input:
		pkls = expand("temp/{{TECH}}_{ID}.pkl", ID=IDs),
	output:
		pkl = "results/{TECH}.pkl"	
	resources:
		mem=32,
	run:
		dfs = []
		total = len(input["pkls"])
		for idx, pkl in enumerate(input["pkls"]):
			dfs.append(pd.read_pickle(pkl) )
			sys.stderr.write("\r{:.02%}".format( (idx+1)/total ) )
		df = pd.concat(dfs, ignore_index=True)
		df.to_pickle(output["pkl"])


rule merge_mat:
	input:
		npys = expand("temp/{{TECH}}_{ID}.mat.npy", ID=IDs),
	output:
		npy = "results/{TECH}.mat.npy"	
	resources:
		mem=32,
	run:
		
		mat = np.load(input["npys"][0])
		out = np.zeros(mat.shape)
		for mat_f in input["npys"]:
			mat = np.load(mat_f)
			out += mat
		np.save(output["npy"], out)
	






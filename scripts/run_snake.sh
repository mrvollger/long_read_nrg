#!/usr/bin/env bash

jobNum=1500
waitTime=60 # this really needs to be 60 on our cluster :(
retry=2 # numer of times to retry the pipeline if it failes
# I allow a retry becuase sometimes even the really long waittime is not enough,
# and the files are actaully there

#
# QSUB parameters, these are only the defualts, they can be changed with params.sge_opts
# Allow snakemake to make directories, I think it slows things down when I done with "waitTime"
#
logDir=logs
mkdir -p $logDir

echo "Running snakemake file: "$1

#
# run snakemake
#
snakemake -p \
        --drmaa " -P eichlerlab \
                -q eichler-short.q \
                -l h_rt=8:00:00  \
                -l mfree={resources.mem}G \
				-pe serial 1 \
                -V -cwd \
                -S /bin/bash" \
		 --drmaa-log-dir $logDir \
        --jobs $jobNum \
        --latency-wait $waitTime \
        --restart-times $retry  \
        -s $@

# generate report 
#snakemake -s $snakefile --report racon_report.html
#-e {log.e} -o {log.o} \



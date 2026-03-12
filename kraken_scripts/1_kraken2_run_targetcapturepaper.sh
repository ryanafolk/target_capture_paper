#!/bin/bash

#classify using Kraken2

# Run in a subdirectory called kraken2_SILVA in your main folder

#Text file with sample names
samples=./samples_targetcapturepaper.txt

#output directory for classified reads
op=./intbio_targetcapturepaper
mkdir $op

#output directory for report files. THIS IS IMPORTANT FOR THE NEXT STEP
rep=./intbio_targetcapturepaper_report
mkdir $rep

#directory with reads. If the sample came from a host organism, I would map out the host reads first.
reads=/mnt/Botbot/sequencing_reads_processed/intbio_reads_NovogeneMay2025

#path to kraken database. For my analysis, I used Silva release 138
db=./16S_SILVA138_k2db

#loop through samples to classify reads
while read -r line; do
	kraken2 --paired --threads 16 --report $rep/${line}.k2.report --db $db $reads/${line}_R1.fq.gz $reads/${line}_R2.fq.gz > $op/${line}_classified.kraken2
	done<$samples
	
# nproc to determine the number of threads available
# also execute top to make sure the system isn't too busy

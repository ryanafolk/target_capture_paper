#!/bin/bash

# Obtain abundances using Bracken

# Run in a subdirectory called kraken2 in your main folder

#Text file with sample names
samples=./samples_targetcapturepaper.txt

#directory with all kraken2 reports for samples.
reports=./intbio_targetcapturepaper_report

#make bracken output directory
op=./intbio_targetcapturepaper_bracken
mkdir $op

#kraken database used to assign taxonomy
db=./16S_SILVA138_k2db

#loop through to make bracken report at all taxonomic levels (up to class). Also disregarding any taxa with less than ten reads mapped. -r is length of reads
while read -r line; do
	bracken -d $db -i $reports/$line.k2.report -o $op/${line}_braken_species -r 150 -l S -t 10; 
	bracken -d $db -i $reports/$line.k2.report -o $op/${line}_braken_genus -r 150 -l G -t 10; 
	bracken -d $db -i $reports/$line.k2.report -o $op/${line}_braken_family -r 150 -l F -t 10; 
	bracken -d $db -i $reports/$line.k2.report -o $op/${line}_braken_order -r 150 -l O -t 10; 
	bracken -d $db -i $reports/$line.k2.report -o $op/${line}_braken_class -r 150 -l C -t 10; 
	done<$samples


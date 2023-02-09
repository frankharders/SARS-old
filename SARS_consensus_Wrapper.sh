#!/bin/bash

## "source activate POPPUNK";

##  Written by Frank Harders adjusted the template of Brian Bushnell
##  Last modified April 26 2021
##  Description: Outer wrapper for processing a bunch of Covid samples in a directory.
##  This is just a template that needs to be modified before use, according to the input file names.

#echo "This template must be modifed before use according to your file names."
#exit

##  Optionally, get rid of old files first, if the pipeline is being rerun.
rm *.sam.gz *.bam *.bai *.txt *_genome.fa *_adapters.fa *.vcf *.vcf.gz;

## create file containing a ist of sample files.

mkdir -p "$PWD"/RAWREADS;

## rename original MiSeq files into a proper format 

rename 's/_S[0-9]*//g' *.gz;
rename 's/_001././g' *.gz;


## create a list of sample names for processing 
ls *R1* | cut -f1 -d'_' > samples.txt;

## check whether sample files are correct
cat samples.txt;


## move all files to a directory 
mv *.gz "$PWD"/RAWREADS;

cd "$PWD"/RAWREADS;
ls *R1* | cut -f1 -d'_' > ./../samples.txt;

cd ..;

sh "$PWD/"create_consensus.sh;

sh "$PWD"/makeSummary.sh 1>makeSummary.o 2>&1;




exit 1

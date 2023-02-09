#!/bin/bash
##  Written by Frank Harders adjusted the template of Brian Bushnell
##  Last modified April 26, 2021
##  Description: Summarizes all the runs of "create_consensus.sh" in this directory.
##  Then tars them if you want to send them somewhere. Should be run only after all samples are processed individually.

## Set minimum coverage for variant calls.
MINCOV=5;

##  Specify the viral reference file.
##  NC_045512.fasta contains the SARS-CoV-2 genome, equivalent to bbmap/resources/Covid19_ref.fa
REFPATH="/home/harde004/miniconda3/envs/POPPUNK/opt/bbmap-38.90-1/resources";
REF="$REFPATH"/"NC_045512.fasta";

cp "$PWD"/02_polished/*/*_genome.fa  "$PWD"/genomes;

mkdir -p "$PWD"/sumAnalysis;

cp "$PWD"/02_polished/*/*_deduped_trimclip.sam.gz  "$PWD"/sumAnalysis;
cp "$PWD"/02_polished/*/*basecov_border5.txt  "$PWD"/sumAnalysis;

##  Call variants in multisample mode to solve things. This reports the genotype of *all* samples at any position at which a variant is called in *any* sample.
callvariants.sh "$PWD"/sumAnalysis/*_deduped_trimclip.sam.gz ref="$REF" multisample out=allVars.vcf ow -Xmx4g usebias=f strandedcov minstrandratio=0 maf=0.6 minreads="$MINCOV" mincov="$MINCOV" minedistmax=30 minedist=16 flagnearby > all.samples.variants.log 2>&1;
##  Make a summary of coverage at varous depth cutoffs for all libraries.
summarizecoverage.sh "$PWD"/sumAnalysis/*basecov_border5.txt out=coverageSummary.txt > "$PWD1"/all.samples.summary.log 2>&1;
mkdir -p output;
cp *.sh output;
cp "$PWD"/02_polished/*/*.bam* output;
cp "$PWD"/02_polished/*/*.txt output;
cp "$PWD"/02_polished/*/*.vcf output;
cp "$PWD"/02_polished/*/*genome*.fa output;
rm results.tar;
tar -cf results.tar output;


exit 1

#!/bin/bash

source activate POPPUNK

##  Written by Frank Harders adjusted the template of Brian Bushnell
##  Last modified April 26 2021
##  Description:  Calls SARS-CoV-2 variants from Illumina amplicon data.
##                This script assumes input data is paired-end.

## Create nescessary directories

workdir="$PWD";

mkdir -p "$workdir"/02_polished "$workdir"/LOGS "$workdir"/TEMP "$workdir"/dbMAPPING "$workdir"/genomes;


##  Set minimum coverage for genotype calls.
##  Areas below this depth will be set to N in the consensus genome.
MINCOV=5;

##  Specify the viral reference file.
##  NC_045512.fasta contains the SARS-CoV-2 genome, equivalent to bbmap/resources/Covid19_ref.fa
REFPATH="/home/harde004/miniconda3/envs/POPPUNK/opt/bbmap-38.90-1/resources";
REF="$REFPATH"/"NC_045512.fa";

echo "$REF";

## PCR primers
## Specify the location of the pimerset used for amplification
## path of used primer location is linked to the environment used for analysis (#LINE 3)
## /home/harde004/miniconda3/envs/POPPUNK/pkgs/bbmap-38.90-h1296035_0/opt/bbmap-38.90-0/resources
PRIMERPATH="/home/harde004/miniconda3/envs/POPPUNK/pkgs/bbmap-38.90-h1296035_0/opt/bbmap-38.90-0/resources";
PRIMERSa="$PRIMERPATH"/"artic3.fasta";
PRIMERSq="$PRIMERPATH"/"qiaseq.fasta";
PRIMERSn="$PRIMERPATH"/"nimagen.fasta";


#MINEDISTMAX=30;
#MINEDIST=16
#MINALLELFRAC=0.05;
#MINQUALMAX=15;
#MINSCORE=15;
#PLOIDY=10;

count0=1;
countS=$(cat samples.txt | wc -l);

#ls ./references/*.fa > used_refs.txt;

#RefsIn=used_refs.txt;

#count1=1;
#countR=$(cat $RefsIn | wc -l);

#REF=./references/MT019532.fa;

#while [ $count1 -le $countR ];do

#	REF=$(cat $RefsIn | awk 'NR=='$count1 );

#mkdir -p ./dbMAPPING/"$REF"/;
#OUT=./dbMAPPING/"$REF"/;

if [ -f "$REF" ] && [ -f "$PRIMERSa" ] ;then 

echo -e "reference file is present";
echo -e "primer file is present";




		while [ $count0 -le $countS ];do

		SAMPLE=$(cat samples.txt | awk 'NR=='$count0 );

		LOG1=./LOGS/1-"$SAMPLE".reformat.log; # done
		LOG2=./LOGS/2-"$SAMPLE".viralreads.log; # done
		LOG3=./LOGS/3-"$SAMPLE".qualCalibration.log; # done
		LOG4=./LOGS/4-"$SAMPLE".adapterPresent.log; # done
		LOG5=./LOGS/5-"$SAMPLE".clumpify.log; # done		
		LOG6=./LOGS/6-"$SAMPLE".adapter-trimming.log; # done	
		LOG7=./LOGS/7-"$SAMPLE".QpcrPrimer-trimming.log; # done
		LOG8=./LOGS/8-"$SAMPLE".alignReference.log; # done
		LOG9=./LOGS/9-"$SAMPLE".dedupe.log; # done
		LOG10=./LOGS/10-"$SAMPLE".deleteChimerisReads.log; # done
		LOG11=./LOGS/11-"$SAMPLE".delete.junkReads.log; # done
		LOG12=./LOGS/12-"$SAMPLE".trimSoftClip.log; # done
		LOG13=./LOGS/13-"$SAMPLE".callvariantsSam.log; # done
		LOG14=./LOGS/14-"$SAMPLE".pileup.log; # done	
		LOG15=./LOGS/15-"$SAMPLE".applyvariants.log; # done		
		LOG16=./LOGS/16-"$SAMPLE".callvariantsConsensus.log; # done


mkdir -p "$workdir"/02_polished/"$SAMPLE";
Filein="$workdir"/02_polished/"$SAMPLE"/"$SAMPLE";

##  This line is in case the script is being re-run, to clear the old output.
rm "$Filein"*.sam.gz "$Filein"*.bam "$Filein"*.bai "$Filein"*.txt "$Filein"*.fa "$Filein"*.vcf;


		echo -e "$SAMPLE";
		echo -e "used reference $REF";

		# 1
		##  If data is paired in twin files, interleave it into a single file.
		##  Otherwise, skip this step.
		##  In this case, the files are assumed to be named "Sample1_R1.fq.gz" and "Sample1_R2.fq.gz"
		reformat.sh in="$Filein"_R#.fq.gz out="$Filein".fq.gz ow=t > "$LOG1" 2>&1;

		# 2
        ##  Split into Covid and non-Covid reads if this has not already been done.
        ##  This step can be skipped if non-Covid was already removed.
		bbduk.sh ow -Xmx1g in="$Filein".fq.gz ref="$REF" outm="$Filein"_viral.fq.gz outu="$Filein"_nonviral.fq.gz k=25> "$LOG2" 2>&1;

		# 3
		##  Recalibrate quality scores prior to any trimming.
		##  *** Requires a recalibration matrix in the working directory (see recal.sh for details). ***
		##  This step is optional but useful for Illumina binned quality scores.
		bbduk.sh in="$Filein"_viral.fq.gz out="$Filein"_recal.fq.gz recalibrate -Xmx1g ow > "$LOG3" 2>&1;

		# 4
		##  Discover adapter sequence for this library based on read overlap.
		##  You can examine the adapters output file afterward if desired;
		##  If there were too few short-insert pairs this step will fail (and you can just use the default Illumina adapters).
		bbmerge.sh in="$Filein"_recal.fq.gz outa="$Filein"_adapters.fa ow reads=1m > "$LOG4" 2>&1;

		# 5
		##  Remove duplicates by sequence similarity.
		##  This is more memory-efficient than dedupebymapping.
		clumpify.sh in="$Filein"_recal.fq.gz out="$Filein"_clumped.fq.gz zl=9 dedupe s=2 passes=4 -Xmx31g > "$LOG5" 2>&1;

		# 6
		##  Perform adapter-trimming on the reads.
		##  Also do quality trimming and filtering.
		##  If desired, also do primer-trimming here by adding, e.g., 'ftl=20' to to trim the leftmost 20 bases.
		##  If the prior adapter-detection step failed, use "ref=adapters"
		bbduk.sh in="$Filein"_clumped.fq.gz out="$Filein"_trimmed.fq.gz minlen=60 ktrim=r k=21 mink=9 hdist=2 hdist2=1 ref="$Filein"_adapters.fa altref=adapters maq=14 qtrim=r trimq=10 maxns=0 tbo tpe ow -Xmx1g ftm=5 > "$LOG6" 2>&1;

		# 7
		##  Trim artic3 primers, if present.
		##  Disable this line if you are not amplifying with primers.
		bbduk.sh in="$Filein"_trimmed.fq.gz out="$Filein"_trimmed2.fq.gz ref="$PRIMERSa" ktrim=l restrictleft=30 k=22 hdist=3 qhdist=1 rcomp=f mm=f > "$LOG7" 2>&1;

		# 8
		##  Align reads to the reference.
		##  Local flag is due to primer-amplification-induced anomalies at read ends;
		##  for randomly-sheared, unamplified data, "local" should be omitted.
		bbmap.sh ref="$REF" in="$Filein"_trimmed2.fq.gz outm="$Filein"_mapped.sam.gz nodisk local maxindel=500 -Xmx4g ow k=12 > "$LOG8" 2>&1;

		# 9
		##  Deduplicate based on mapping coordinates.
		##  Note that if you use single-ended amplicon data, you will lose most of your data here.
		dedupebymapping.sh in="$Filein"_mapped.sam.gz out="$Filein"_deduped.sam.gz -Xmx31g ow > "$LOG9" 2>&1;

		# 10
		##  Remove junk reads with unsupported unique deletions; these are often chimeric.
		filtersam.sh ref="$REF" ow in="$Filein"_deduped.sam.gz out="$Filein"_filtered.sam.gz mbad=1 del sub=f mbv=0 -Xmx4g > "$LOG10" 2>&1;

		# 11
		##  Remove junk reads with multiple unsupported unique substitutions; these are often junk, particularly on Novaseq.
		##  This step is not essential but reduces noise.
		filtersam.sh ref="$REF" ow in="$Filein"_filtered.sam.gz out="$Filein"_filtered2.sam.gz mbad=1 sub mbv=2 -Xmx4g > "$LOG11" 2>&1;

		# 12
		##  Trim soft-clipped bases.
		bbduk.sh in="$Filein"_filtered2.sam.gz trimclip out="$Filein"_trimclip.sam.gz -Xmx1g ow > "$LOG12" 2>&1;

		# 13
		##  Call variants from the sam files.
		##  The usebias=f/minstrandratio=0 flags are necessary due to amplicon-induced strand bias,
		##  and should be removed if the data is exclusively shotgun/metagenomic or otherwise randomly fragmented.
		callvariants.sh in="$Filein"_trimclip.sam.gz ref="$REF" out="$Filein"_vars.vcf -Xmx4g ow strandedcov usebias=f minstrandratio=0 maf=0.6 minreads="$MINCOV" mincov="$MINCOV" minedistmax=30 minedist=16 flagnearby  > "$LOG13" 2>&1;

		# 14
		##  Calculate reduced coverage as per CallVariants defaults (ignoring outermost 5bp of reads).
		pileup.sh in="$Filein"_trimclip.sam.gz basecov="$Filein"_basecov_border5.txt -Xmx4g ow border=5  > "$LOG14" 2>&1;

		# 15
		##  Generate a mutant reference by applying the detected variants to the reference.
		##  This is essentially the reference-guided assembly of the strain.
		##  Also changes anything below depth MINCOV to N (via the mindepth flag).
		##  Does not apply indels below MINCOV
		applyvariants.sh in="$REF" out="$Filein"_genome.fa vcf="$NAME"_vars.vcf basecov="$Filein"_basecov_border5.txt ow mindepth="$MINCOV" > "$LOG14" 2>&1;

		# 16
		##  Call variants from the sam files.
		##  The usebias=f/minstrandratio=0 flags are necessary due to amplicon-induced strand bias,
		##  and should be removed if the data is exclusively shotgun/metagenomic or otherwise randomly fragmented.
		callvariants.sh in="$Filein"_trimclip.sam.gz ref="$Filein"_genome.fa out="$Filein"_vars.vcf -Xmx4g ow strandedcov usebias=f minstrandratio=0 maf=0.05 minreads="$MINCOV" mincov="$MINCOV" minedistmax=30 minedist=16 flagnearby  > "$LOG16" 2>&1;

		# 17
		##  Make bam/bai files; requires samtools to be installed.
		##  This step is only necessary for visualization, not variant-calling.
		samtools view -bShu "$Filein"_trimclip.sam.gz | samtools sort -m 2G -@ 3 - -o "$Filein"_sorted.bam;
		samtools index "$Filein"_sorted.bam;


######
#
#		clumpify.sh in=./02_polished/"$SAMPLE".fq.gz out=./02_polished/"$SAMPLE"_clumped.fq.gz zl=9 dedupe s=2 passes=4 -Xmx31g > "$LOG2" 2>&1;
#		echo 'bbduk_1';
#		# 3
#		bbduk.sh in=./02_polished/"$SAMPLE"_clumped.fq.gz out=./02_polished/"$SAMPLE"_trimmed.fq.gz minlen=60 ktrim=r k=21 mink=9 hdist=2 hdist2=1 ref=./adapter.fa altref=adapters maq=14 qtrim=r trimq=10 maxns=0 tbo tpe ow -Xmx1g ftm=5 > "$LOG3" 2>&1;
#		echo 'bbmap'
#		# 4
#		bbmap.sh ref="$REF" in=./02_polished/"$SAMPLE"_trimmed.fq.gz outm="$OUT"/"$SAMPLE"_mapped.sam.gz nodisk local maxindel=500 -Xmx4g ow k=12 > "$LOG4" 2>&1;
#		echo 'dedupe'
#		# 5
#		dedupebymapping.sh in="$OUT"/"$SAMPLE"_mapped.sam.gz out="$OUT"/"$SAMPLE"_deduped.sam.gz -Xmx31g ow > "$LOG5" 2>&1;
#		echo 'filtersam_1'
#		# 6
#		filtersam.sh ref="$REF" ow in="$OUT"/"$SAMPLE"_deduped.sam.gz out="$OUT"/"$SAMPLE"_filtered.sam.gz mbad=1 del sub=f mbv=0 -Xmx4g > "$LOG6" 2>&1;
#		echo 'filtersam_2'
#		# 7
#		filtersam.sh ref="$REF" ow in="$OUT"/"$SAMPLE"_filtered.sam.gz out="$OUT"/"$SAMPLE"_filtered2.sam.gz mbad=1 sub mbv=2 -Xmx4g > "$LOG7" 2>&1;
#		echo 'bbduk_2'
#		# 8
#		bbduk.sh in="$OUT"/"$SAMPLE"_filtered2.sam.gz trimclip out="$OUT"/"$SAMPLE"_trimclip.sam.gz -Xmx1g ow > "$LOG8" 2>&1;
#		echo 'callvariants'
#		# 9
#		callvariants.sh in="$OUT"/"$SAMPLE"_trimclip.sam.gz ref="$REF" out="$OUT"/"$SAMPLE"_vars.vcf -Xmx4g ow strandedcov usebias=f minstrandratio=0 maf=0.05 minreads="$MINCOV" mincov="$MINCOV" minedistmax="$MINEDISTMAX" minedist="$MINEDIST" flagnearby > "$LOG9" 2>&1;
#		echo 'pileup'
#		# 10
#		pileup.sh in="$OUT"/"$SAMPLE"_trimclip.sam.gz basecov="$OUT"/"$SAMPLE"_basecov_border5.txt -Xmx4g ow border=5 > "$LOG10" 2>&1;
#		echo 'applyvariants'
#		# 11
#		applyvariants.sh in="$REF" out="$OUT"/"$SAMPLE"_genome.fa vcf="$OUT"/"$SAMPLE"_vars.vcf basecov="$OUT"/"$SAMPLE"_basecov_border5.txt ow mindepth="$MINCOV" > "$LOG11" 2>&1;
#		echo 'mapping original reads to "consensus sequence"'
#		# 12
#		bbrename.sh in="$OUT"/"$SAMPLE"_genome.fa out="$OUT"/"$SAMPLE"_genome.consensus.fa prefix="$SAMPLE";
#
#		REF1="$OUT"/"$SAMPLE"_genome.consensus.fa;
#		#13
#		bbmap.sh ref="$REF1" in=./02_polished/"$SAMPLE"_trimmed.fq.gz outm="$OUT"/"$SAMPLE"_mapped.consensus.sam.gz nodisk local maxindel=500 -Xmx4g ow k=12 > "$LOG12" 2>&1;
#		#14
#		samtools view -bShu "$OUT"/"$SAMPLE"_mapped.consensus.sam.gz | samtools sort -m 2G -@ 3 - -o "$OUT"/"$SAMPLE"_sorted.consensus.bam;
#		samtools index "$OUT"/"$SAMPLE"_sorted.consensus.bam;
#		#15
#		echo 'call variants on consensus sequence';
#		callvariants.sh in="$OUT"/"$SAMPLE"_sorted.consensus.bam out="$OUT"/"$SAMPLE".genome.consensus.vcf ref="$REF1" minallelefraction="$MINALLELFRAC" minqualitymax="$MINQUALMAX" minscore="$MINSCORE" ploidy="$PLOIDY" ow strandedcov usebias=f minstrandratio=0 maf=0.05 minreads="$MINCOV" mincov="$MINCOV" minedistmax="$MINEDISTMAX" minedist="$MINEDIST" flagnearby > "$LOG13" 2>&1;

		count0=$((count0+1));

	done

#count1=$((count1+1));


else 

echo "reference or primer file is missing from the original location";
echo "check if the proper environment is selected for analysis";
echo "make sure that these files are present at this location $REFPATH"/"NC_045512.fa or $PRIMERPATH"/"artic3.fasta";


fi 








#done 
exit 1


#Variant-Calling Cutoffs:
#minreads=2              (minad) Ignore variants seen in fewer reads.
#maxreads=BIG            (maxad) Ignore variants seen in more reads.
#mincov=0                Ignore variants in lower-coverage locations.
#maxcov=BIG              Ignore variants in higher-coverage locations.
#minqualitymax=15        Ignore variants with lower max base quality.
#minedistmax=20          Ignore variants with lower max distance from read ends.
#minmapqmax=0            Ignore variants with lower max mapq.
#minidmax=0              Ignore variants with lower max read identity.
#minpairingrate=0.1      Ignore variants with lower pairing rate.
#minstrandratio=0.1      Ignore variants with lower plus/minus strand ratio.
#minquality=12.0         Ignore variants with lower average base quality.
#minedist=10.0           Ignore variants with lower average distance from ends.
#minavgmapq=0.0          Ignore variants with lower average mapq.
#minallelefraction=0.1   Ignore variants with lower allele fraction.  This
#                        should be adjusted for high ploidies.
#minid=0                 Ignore variants with lower average read identity.
#minscore=20.0           Ignore variants with lower Phred-scaled score.
#clearfilters            Clear all filters.  Filter flags placed after
#                        the clearfilters flag will still be applied.







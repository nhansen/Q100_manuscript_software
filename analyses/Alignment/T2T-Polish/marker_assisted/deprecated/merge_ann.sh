#!/usr/bin/env bash

######################################################################
#  PUBLIC DOMAIN NOTICE
#
#  This software is "United States Government Work" under the terms of the United
#  States Copyright Act. It was written as part of the authors' official duties
#  for the United States Government and thus cannot be copyrighted. This software
#  is freely available to the public for use without a copyright
#  notice. Restrictions cannot be placed on its present or future use.
#
#  Although all reasonable efforts have been taken to ensure the accuracy and
#  reliability of the software and associated data, the National Human Genome
#  Research Institute (NHGRI), National Institutes of Health (NIH) and the
#  U.S. Government do not and cannot warrant the performance or results that may
#  be obtained by using this software or data. NHGRI, NIH and the U.S. Government
#  disclaim all warranties as to performance, merchantability or fitness for any
#  particular purpose.
#
#  Please cite the authors in any work or product based on this material.
######################################################################

if [[ "$#" -lt 4 ]]; then
	echo "Usage: single_copy_filter.sh alignment target asm marker.meryl [len_filt]"
	echo
	echo "Filters alignment based on single-copy kmers"
	echo -e "\talignment: input bam or cram file, containing *only* target region (chr)."
	echo -e "\ttarget: target region (chr) to process"
	echo -e "\tasm: target fasta file"
	echo -e "\tmarker.meryl: meryl db containing marker kmers. This scripts generates chr specific markers."
	echo -e "\tlen_filt: length filter in kb. Read from LEN_FILT. OPTIONAL"
	exit -1
fi


## Dependency: samtools, samToAlignment, meryl v1.0, subsetSamByKmers.py, samToErrorRate, SubFile, IGVTools

alignment=$1    # input bam / cram file 
target=$2	# chrX
asm=$3		# asm.fasta
single=$4	# single-copy.meryl
len_filt=$5	# length filter

cores=$SLURM_CPUS_PER_TASK
if [ x$cores == "x" ]; then
	echo "Use 16 cores by default"
	cores=16
fi

if [ -z $len_filt ]; then
	len_filt=`cat LEN_FILT`
fi
len_filt=$((len_filt*1000))


SCRIPT=`cat SCRIPT`
module load samtools/1.9

set -e

if [ ! -e $asm.fai ]; then
	echo "init.sh not executed"
	exit 0
fi
if [ ! -e $target.fa ]; then
	echo "init.sh not executed"
	exit 0
fi

if [ -e $target.markersandlength.cram ]; then
	echo "Already done"
	exit 0
else
	if [ ! -s $target.aligned.fasta ]; then
		#Following steps are run after job arrays (main.sh) have finished

		echo "Merge to $target.aligned.fasta"
		cat split.$target.srt_id.*.aligned.fasta > $target.aligned.fasta
		echo
	fi

	if [ ! -s $target.alignment.posCount ]; then
		echo "Prepare $target.single.meryl"
		/data/Phillippy/tools/meryl-tip/build/bin/meryl count k=21 $target.fa output $target.meryl
		/data/Phillippy/tools/meryl-tip/build/bin/meryl intersect $target.meryl $single output $target.single.meryl
		/data/Phillippy/tools/meryl-tip/build/bin/meryl-lookup -existence -memory 12 -sequence $target.aligned.fasta -mers $target.single.meryl | awk '$NF>0 {print $1"\t"$NF}' > $target.alignment.posCount
		echo
	fi

	echo "Subset by markers"
	if [ ! -s $target.markers.cram ]; then
		samtools view -h -O sam -@$cores $alignment > $target.sam
		echo "python $SCRIPT/src/subsetSamByKmers.py $target.alignment.posCount $target.sam > $target.markers.sam"
		python $SCRIPT/src/subsetSamByKmers.py $target.alignment.posCount $target.sam > $target.markers.sam
		samtools view -@$cores -O cram -o $target.markers.cram --reference=$asm $target.markers.sam
		$tools/IGVTools/igvtools count $target.markers.cram $target.markers.tdf $asm.fai
	fi
	echo

	echo "# hard filter alignments < $len_filt kb to $target.filtered.sam"
	if [ -s $target.filtered.sam ]; then
		echo "*** Found $target.filtered.sam. Skipping this step. ***"
	else
		cat $target.header > $target.filtered.sam
		samtools view -@$cores $target.markers.sam |awk -v l=$len_filt '{if (length($10) >= l) print $0}' >> $target.filtered.sam
	fi
	echo

	echo "# also filter by alignment legnth and 75% idy"
	if [ -s $target.filteredList ]; then
		echo "*** Found $target.filteredList. Skipping this step. ***"
	else
		$SCRIPT/src/samToErrorRate $target.filtered.sam $asm \
		| awk '{if ($9 == 0) print $0; else print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t1\t"$12-$11"\t"$12-$10"\t"$12"\t"$13"\t"$14"\t"$15"\t"$16"\t"$17}' \
		| awk -v l=$len_filt '{if ($3 <= -1*l && $4 >= 75) print $target}'| awk -v rid="" '{if (rid!=$1) {print $1; rid=$1}}' > $target.filteredList
	fi
	cat $target.header > $target.tmp.sam
	samtools view -@$cores $target.filtered.sam > $target.tmp2.sam
	java -cp $SCRIPT/src/ SubFile $target.filteredList $target.tmp2.sam >> $target.tmp.sam
	mv $target.tmp.sam $target.filtered.sam
	rm $target.tmp2.sam
fi

echo "
# generate $target.markersandlength.cram"
samtools view -@${cores} -O cram -o $target.markersandlength.cram --reference=$asm $target.filtered.sam
$tools/IGVTools/igvtools count $target.markersandlength.cram $target.markersandlength.tdf $asm.fai

echo "
# Index"
samtools index $target.markersandlength.cram

echo "
# Cleanup"
rm $target.filtered.sam $target.sam $target.markers.sam $target.bam split.$target.srt_id.* $target.align* $target.filteredList

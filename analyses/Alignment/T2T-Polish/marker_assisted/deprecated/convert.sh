#!/usr/bin/env bash

#SBATCH --gres=lscratch:5

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

if [[ "$#" -lt 2 ]]; then
	echo "Usage: convert.sh target asm [arr_id]"
	echo
	echo "Convert alignments to position preserved .fasta file"
	echo -e "\ttarget: target (chr) to process"
	echo -e "\tasm: assembly fasta file"
	echo -e "\t[arr_id]: array id to process. DEFAULT=\$SLURM_ARRAY_TASK_ID"
	exit -1
fi


## Dependency: samToAlignments

target=$1	# chrX
asm=$2		# asm.fasta
i=$SLURM_ARRAY_TASK_ID

cores=$SLURM_CPUS_PER_TASK
if [ x$cores == "x" ]; then
	echo "Use 16 cores by default"
	cores=16
fi

if [[ -z $i ]]; then
	if [[ -z $3 ]]; then
		echo "No job array id provided. Submit through a job array or set [arr_id]"
		exit -1
	fi
	i=$3
fi

SCRIPT=`cat SCRIPT`

set -e

if [ ! -e $asm.fai ]; then
  echo "init.sh not executed"
  exit 0
fi

if [ -e $target.markersandlength.cram ]; then
  echo "Already done"
  exit 0
elif [ ! -s split.$target.aligned.fasta ]; then
  split=split.$target.srt_id.$i
  echo "Processing file $split"
  cat $target.header > $split.tmp.sam
  # Sort alignments by flag for each read ID" ## commenting out as printing these logs take much longer
  java -jar -Xmx12g $SCRIPT/src/samSortByFlag.jar $split >> $split.tmp.sam
  echo -e "\nsamToAlignment $split.tmp.sam $asm"
  $SCRIPT/src/samToAlignment $split.tmp.sam $asm 2> $split.tmp.err \
    | awk '{if ($9 == 0) { print ">"$1"_"$5"_"$9"_"$10"_"$(NF-1); print $NF } else { print ">"$1"_"$5"_"$9"_"$12-$11"_"$(NF-1); print $NF}}' \
    > $split.aligned.fasta

  echo
  echo "Clean up $split $(ls $split.tmp.*)"
  rm $split $split.tmp.*

  echo "Output file generated as $split.aligned.fasta"
else
  echo "split.$target.aligned.fasta already exists. Nothing to do."
fi

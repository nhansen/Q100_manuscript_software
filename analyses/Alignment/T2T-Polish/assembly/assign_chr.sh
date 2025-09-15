#!/bin/sh

if [[ "$#" -lt 2 ]]; then
  echo "Usage: ./assign_chr.sh asm.fa ref.fa"
  echo "Assign chromosome names in asm.fa based on ref.fa"
  echo "Output will be written to assembly.refOriented.haplotype.fa"
  exit -1
fi

asm=$1
ref=$2

module load samtools/1.21
module load mashmap/3.1.1

set -x

if [[ ! -s $asm.fai ]]; then
  samtools faidx -@12 $asm
fi

if [[ ! -s assembly-ref.norm.mashmap ]]; then
  mashmap \
    --ref $ref \
    --query $asm \
    --perc_identity 95 \
    --segLength 100000 \
    --threads 12 \
    --output assembly-ref.norm.mashmap
fi

if grep -q mat- $asm ; then
    label1="pat-"
    label2="mat-"
elif grep -q h1tg $asm ; then
    label1="h1tg"
    label2="h2tg"
elif grep -q contig- $asm ; then
    label1="contig-"
    label2="none-ignore"
elif grep -q "#1#" $asm ; then
   label1="#1#"
   label2="#2#"
else
    label1="haplotype1-"
    label2="haplotype2-"
fi
echo "label1: $label1 | label2: $label2"

if [ ! -e translation_hap1 -o ! -e translation_hap2 ] ; then
    cat assembly-ref.norm.mashmap | grep "$label1" | \
	  awk '{if ($NF > 99 && $4-$3 > 1000000) print $1"\t"$6"\t"$2"\t"$7}' |\
	    sort | uniq > translation_hap1
    cat assembly-ref.norm.mashmap | grep "$label2" | \
	  awk '{if ($NF > 99 && $4-$3 > 1000000) print $1"\t"$6"\t"$2"\t"$7}' |\
	    sort | uniq > translation_hap2
fi

if [ ! -e chr_completeness_max_hap1 -o ! -e chr_completeness_max_hap2 ] ; then
    cat translation_hap1 | sort -k2,2 | awk '{if ($3 > 15000000) print $0}' |\
	  awk -v LAST="" -v S="" '{if (LAST != $2) { \
	    if (S > 0) print LAST"\t"C"\t"SUM/S*100"\t"MAX/S*100"\t"TIG; SUM=0; MAX=0; C=0; } \
	    LAST=$2; S=$NF; SUM+=$3; if (MAX < $3) MAX=$3; C+=1; TIG=$1} \
	    END {print LAST"\t"C"\t"SUM/S*100"\t"MAX/S*100"\t"TIG;}' |\
		 awk '{print $1"\t"$4}' | sort -nk1,1 -s > chr_completeness_max_hap1
    cat translation_hap2 | sort -k2,2 | awk '{if ($3 > 15000000) print $0}' |\
	  awk -v LAST="" -v S="" '{if (LAST != $2) { \
	    if (S > 0) print LAST"\t"C"\t"SUM/S*100"\t"MAX/S*100"\t"TIG; SUM=0; MAX=0; C=0; } \
	    LAST=$2; S=$NF; SUM+=$3; if (MAX < $3) MAX=$3; C+=1; TIG=$1} \
		END {print LAST"\t"C"\t"SUM/S*100"\t"MAX/S*100"\t"TIG;}' |\
		 awk '{print $1"\t"$4}' | sort -nk1,1 -s > chr_completeness_max_hap2
fi

# this is chromosome assignment
if [[ ! -s assembly.refOriented.fasta ]]; then
   sh $tools/T2T-Polish/assembly/reorientByRef.sh assembly-ref.norm.mashmap > assembly-ref.reorient.tsv
   cat assembly-ref.reorient.tsv | awk '{print $1}' > tmp
   grep -w -v -f tmp $asm.fai | awk '{print $1"\t0\t"$2}' >>  assembly-ref.reorient.tsv
   rm ./tmp

   java -cp /data/korens/devel/utils:. SubFasta assembly-ref.reorient.tsv $asm >  assembly.refOriented.fasta
fi

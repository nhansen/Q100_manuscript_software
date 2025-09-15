#!/bin/sh

module load bedtools/2.31.1

MIN_LEN=200000
if [[ "$#" -ge 2 ]]; then
  MIN_LEN=$2
fi
MIN_FRAC=10
if [[ "$#" -ge 3 ]]; then
   MIN_FRAC=$3
fi
MAPPED=""

IGNORED_BED="$tools/T2T-Polish/assembly/chm13_ignore.bed"

for i in `cat $1 |awk -v M=$MIN_LEN '{if ($2 >= M) print $1}'|sort |uniq`; do
   isMAPPED=`cat $1 |grep -w $i |sed s/id:f://g |awk '{if ($(NF-1) > 0.99) print $6"\t"$8"\t"$9"\t"$0}' | \
     bedtools intersect -f 1 -v -a - -b $IGNORED_BED  | wc -l`
   if [ $isMAPPED -eq 0 ]; then
      echo "WARNING: no mapping for $i" >&2
	  continue
   fi
   A=`cat $1 |grep -w $i |sed s/id:f://g |awk '{if ($(NF-1) > 0.99) { print $6"\t"$8"\t"$9"\t"$0; } }' | \
     bedtools intersect -f 1 -v -a - -b $IGNORED_BED | \
       awk '{print $4"\t"$5"\t"$7-$6"\t"$8"\t"$9"\t"$12-$11}' | \
         sort -sk5,5 |\
           awk -v PREV="" '{if (PREV != $5) {if (PREV != "") print PREV"\t"$1"\t"FWD/LEN*100"\t"REV/LEN*100"\t"LEN; PREV=$5; LEN=$2; FWD=0; REV=0; } if ($4 == "-") REV+=$3;  else FWD+=$3;  } END { print PREV"\t"$1"\t"FWD/LEN*100"\t"REV/LEN*100"\t"LEN }' | \
             awk -v MF=$MIN_FRAC '{if ($3 > MF && $3 > $4) print $2"\t0\t"$NF"\t"$1"_"$2"\t"$3; else if ($4 > MF) print $2"\t"$NF"\t0\t"$1"_"$2"\t"$4}'`
   isOK=`cat $1 |grep -w $i |sed s/id:f://g |awk '{if ($(NF-1) > 0.99) { print $6"\t"$8"\t"$9"\t"$0; } }' | \
     bedtools intersect -f 1 -v -a - -b $IGNORED_BED | \
       awk '{print $4"\t"$5"\t"$7-$6"\t"$8"\t"$9"\t"$12-$11}' | \
         sort -sk5,5 |awk -v PREV="" '{if (PREV != $5) {if (PREV != "") print PREV"\t"$1"\t"FWD/LEN*100"\t"REV/LEN*100"\t"LEN; PREV=$5; LEN=$2; FWD=0; REV=0; } if ($4 == "-") REV+=$3;  else FWD+=$3;  } END { print PREV"\t"$1"\t"FWD/LEN*100"\t"REV/LEN*100"\t"LEN }' | \
           awk -v MF=$MIN_FRAC '{if ($3 > MF && $3 > $4) print $2"\t0\t"$NF"\t"$1"_"$2; else if ($4 > MF) print $2"\t"$NF"\t0\t"$1"_"$2}' | \
            grep -v "^$" | wc -l`
   if [ $isOK -eq 0 ]; then
      echo "WARNING: could not assign $i" >&2
   elif [ $isOK -gt 1 ]; then
      echo "ERROR: more than one assignment for $i, it has $A" >&2
	  echo $A | tr '\n' ';' | sed 's/;$/\n/g'
   else
      echo $A
   fi
done

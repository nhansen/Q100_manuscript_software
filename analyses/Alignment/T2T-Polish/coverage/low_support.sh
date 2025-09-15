#!/bin/bash

if ! [[ "$#" -eq 4 ]]; then
  echo "Usage: ./low_support.sh in.paf ver platform pattern"
  echo "  in.paf    bam2paf"
  echo "  ver       assembly version prefix"
  echo "  platform  HiFi, ONT, or Hybrid"
  echo "  pattern   path to pattern folder" 
  exit -1
fi

paf=$1
ver=$2
pattern=$4

# link pattern dir
if ! [[ -d pattern ]]; then
  ln -sf $pattern pattern
fi

asm=pattern/$ver.bed               # chr length
telo=pattern/$ver.telo.bed         # telomere annotation
exclude=pattern/$ver.exclude.bed   # region to exclude
error=pattern/$ver.error.bed       # Merqury asm_only track, merged
pattern=pattern/$ver               # seqrequester microsatellite annotation

pre=${paf/.paf/}

COL_LOW="204,0,0"      # red
COL_CLP="153,153,153"  # grey
COL_HIG="153,102,255"  # purple

LOW=`cat low_cutoff.txt`
HIGH=`cat high_cutoff.txt`

set -e
set -o pipefail

export asset=$tools/asset
module load bedtools

if [[ "$3" == "ONT" ]] || [[ "$3" == "ont" ]]; then
  asset_opt=""
  clip_thresh="15"
elif [[ "$3" == "HiFi" ]] || [[ "$3" == "hifi" ]]; then
  asset_opt="-l 0"
  clip_thresh="10"
elif [[ "$3" == "Hybrid" ]] || [[ "$3" == "hybrid" ]]; then
  asset_opt="-l 0"
  clip_thresh="5"
else
  echo "$3 not recognizable. Exit."
  exit -1
fi

echo "#### $3 Mode. LOW: $LOW, HIGH: $HIGH ####"

echo "# Collect low coverage with Asset; $asset_opt"
set -x

#:<<'END'
$asset/bin/ast_pb $asset_opt -m $LOW -M 10000000 $paf > $pre.bed
rm pb.cov.wig

echo "# Extract low coverage"
bedtools subtract -a $asm -b $pre.bed | bedtools merge -d 5000 -i - | awk -v col=$COL_LOW '{print $0"\tLow\t100\t.\t"$2"\t"$3"\t"col}' > $pre.low.bed

echo "# Collect high coverage"
$asset/bin/ast_pb $asset_opt -m 0 -M $HIGH $paf > $pre.bed
mv pb.cov.wig $pre.cov.asset.wig

echo "# Extract high coverage"
bedtools subtract -a $asm -b $pre.bed | bedtools merge -d 5000 -i - | awk -v col=$COL_HIG '{print $0"\tHigh\t100\t.\t"$2"\t"$3"\t"col}' > $pre.high.bed

echo "
# Collect clipped.bed using clip_*.wig from collect_summary.sh"
java -jar -Xmx1g $tools/T2T-Polish/paf_util/wigToBed.jar $pre.clip_abs.wig  | awk -v t=$clip_thresh '$NF>t' | bedtools merge -i - > $pre.clip_abs.bed
java -jar -Xmx1g $tools/T2T-Polish/paf_util/wigToBed.jar $pre.clip_norm.wig | awk -v t=$clip_thresh '$NF>t' | bedtools merge -i - > $pre.clip_norm.bed
cat $pre.clip_*.bed | bedtools sort -i - | bedtools merge -i - | awk -v col=$COL_CLP '{print $0"\tClipped\t100\t.\t"$2"\t"$3"\t"col}' > clipped.bed

window=10000

for coverage in $pre.low $pre.high
do
  bedtools window -w $window -c -a $coverage.bed    -b $pattern/$ver.microsatellite.GA.128.gt80.bed > $coverage.ga.bed
  bedtools window -w $window -c -a $coverage.ga.bed -b $pattern/$ver.microsatellite.TC.128.gt80.bed > $coverage.tc.bed
  bedtools window -w $window -c -a $coverage.tc.bed -b $pattern/$ver.microsatellite.GC.128.gt80.bed > $coverage.gc.bed
  bedtools window -w $window -c -a $coverage.gc.bed -b $pattern/$ver.microsatellite.AT.128.gt80.bed > $coverage.at.bed
  bedtools intersect -c -a $coverage.at.bed -b $error > $coverage.err.bed
done

coverage=$pre.low
cat $coverage.err.bed | awk -v col_low=$COL_LOW '{ if ($NF>0) print $1"\t"$2"\t"$3"\tLow_Qual\t"$5"\t"$6"\t"$7"\t"$8"\t"col_low; else if ($10+$11 > 0) print $1"\t"$2"\t"$3"\tLow_GA/TC\t"$5"\t"$6"\t"$7"\t"$8"\t153,153,255";  else if ($12 > 0) print $1"\t"$2"\t"$3"\tLow_GC\t"$5"\t"$6"\t"$7"\t"$8"\t204,153,255";  else if ($13 > 0) print $1"\t"$2"\t"$3"\tLow_AT\t"$5"\t"$6"\t"$7"\t"$8"\t204,153,255";  else print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9; }' - > $coverage.pattern.bed

coverage=$pre.high
cat $coverage.err.bed | awk -v col_low=$COL_HIG '{ if ($NF>0) print $1"\t"$2"\t"$3"\tHigh_Qual\t"$5"\t"$6"\t"$7"\t"$8"\t"col_low; else if ($10+$11 > 0) print $1"\t"$2"\t"$3"\tHigh_GA/TC\t"$5"\t"$6"\t"$7"\t"$8"\t153,153,255";  else if ($12 > 0) print $1"\t"$2"\t"$3"\tHigh_GC\t"$5"\t"$6"\t"$7"\t"$8"\t204,153,255";  else if ($13 > 0) print $1"\t"$2"\t"$3"\tHigh_AT\t"$5"\t"$6"\t"$7"\t"$8"\t204,153,255";  else print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9; }' - > $coverage.pattern.bed

cat $pre.low.pattern.bed $pre.high.pattern.bed > $pre.low_high.pattern.bed
coverage=$pre.low_high

echo "
# Collect clipped flanking with $coverage.pattern.bed"
cut -f1-3 $coverage.pattern.bed | bedtools window -w 5 -u -a clipped.bed -b - > clipped_only.bed
cat clipped_only.bed $coverage.pattern.bed | bedtools sort -i - | bedtools subtract -A -a - -b $exclude | bedtools subtract -A -a - -b $telo > $pre.issues.bed

echo "
# Clean up"
rm clipped_only.bed $pre.low_high.pattern.bed
for coverage in $pre.low $pre.high
do
  rm $coverage.ga.bed $coverage.tc.bed $coverage.gc.bed $coverage.at.bed $coverage.err.bed $coverage.pattern.bed
done

echo "
# Collect clippings not reported in issues.bed"
bedtools subtract -A -a clipped.bed -b $pre.issues.bed | bedtools subtract -A -a - -b $exclude | bedtools subtract -A -a - -b $telo > $pre.clipped.bed




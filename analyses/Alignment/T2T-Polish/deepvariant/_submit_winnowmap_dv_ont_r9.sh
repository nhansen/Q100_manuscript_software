#!/bin/bash

if [[ "$#" -lt 3 ]]; then
  echo "Usage: \$tools/T2T-Polish/deepvariant/_submit_winnowmap_dv_ont_r9.sh ref.fa out-prefix mq [-y] [jid]"
  echo "  ref         reference fasta. REQUIRES .fai file in the same place"
  echo "  out-prefix  output prefix"
  echo "  mq          minimum mapping quality requirements. Use 0 for all-to-dip, 1 for all-to-hap alignment."
  echo "  -y          OPTIONAL. set for passing methylation tags."
  echo "  jid         OPTIONAL. Slurm winnowmap jobid to wait, if it is already running or completed."
  echo
  echo "  output      output files will be generated under dv_ONT_R9_MINMQ\${mq}_chr."
  exit -1
fi

ref=$1
out=$2
mq=$3
map="map-ont $4"
jid=$5

if [[ -z $jid ]]; then
  $tools/T2T-Polish/winnowmap/_submit.sh $ref $out $map
  jid=`cat filt.jid | tail -n1`
else
  echo "Skip submitting winnowmap. Wait for $jid"
fi

set -x
$tools/T2T-Polish/deepvariant/_submit_ont_r9_pepper_margin_dv.sh $out.pri.bam $ref $mq $jid >> ont_dv.mrg.jid
set +x

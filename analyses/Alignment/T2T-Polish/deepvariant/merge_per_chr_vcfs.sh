#!/bin/bash

if [[ "$#" -lt 2 ]]; then
  echo "Usage: ./merge_per_chr_vcfs.sh in_dir_prefix out"
  echo "  in_dir_prefix  input directory prefix. something like dv_ONT_R9_"
  echo "  out            output vcf prefix."
  echo
  echo "Output file will be generated as \$out.vcf.gz"
  exit -1
fi

in=$1
out=$2
cpus=$SLURM_CPUS_PER_TASK

ls ${in}*/PEPPER_MARGIN_DEEPVARIANT_FINAL_OUTPUT.vcf.gz > r9_files_to_mrg.list

module load bcftools
module load deepvariant

set -x
bcftools concat -D -a --threads $cpus --no-version -Oz -o $out.vcf.gz -f r9_files_to_mrg.list
bcftools index  $out.vcf.gz
vcf_stats_report --input_vcf $out.vcf.gz --outfile_base $out

#! /bin/bash

echo "Usage: ./step2.sh"

if [[ ! -s MODE ]]; then
  echo "Requires MODE text file to be present. Exit."
  exit -1
fi

# Update to DeepVariant v1.6.1 on Sep. 26 2024
module load deepvariant/1.6.1 || exit 1

ulimit -u 32768
set -o pipefail
set -e

MODE=`cat MODE`
OUT=dv_$MODE
mode=`echo $MODE | awk '{print tolower($0)}'`
MODEL="/opt/models/$mode/model.ckpt"
N_SHARD=`cat N_SHARD` # Must match as in step1
CALL_VARIANTS_OUTPUT="$OUT/call_variants_output.tfrecord.gz"

set -x

call_variants \
  --outfile "${CALL_VARIANTS_OUTPUT}" \
  --examples "$OUT/examples/tfrecord@${N_SHARD}.gz" \
  --checkpoint "${MODEL}" \
  || exit 1

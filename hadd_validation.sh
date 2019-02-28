#/bin/csh

set input=$1
set mode=$2
set period=$3

rm logs/hadd_validation_${input}_mode_${mode}_$period
hadd -f output/${input}.${period}.mode_${mode}_validation.root output_batch/${input}_${mode}/${period}/file/${input}.${period}.mode_${mode}.* > & logs/hadd_validation_${input}_mode_${mode}_$period

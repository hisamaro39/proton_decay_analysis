#/bin/csh

set input=$1
set mode=$2
set type=$3
set period=$4

rm logs/hadd_${input}_mode_${mode}_${type}_$period
hadd -f output/${input}.${period}.mode_${mode}_${type}.root output_batch/${input}_${mode}/${period}/file/${input}.${period}.mode_${mode}_${type}.* > & logs/hadd_${input}_mode_${mode}_${type}_$period

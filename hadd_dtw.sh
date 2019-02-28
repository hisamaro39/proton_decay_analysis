#/bin/csh

set input=$1
set mode=$2
set period=$3
set dtw=$4

rm logs/hadd_${input}_mode_${mode}_$period
hadd -f output/${input}.${period}.mode_${mode}_distance_to_wall_${dtw}.root output_batch/${input}_${mode}/${period}/file/${input}.${period}.mode_${mode}.* > & logs/hadd_${input}_mode_${mode}_$period

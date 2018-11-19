#/bin/csh

set input=$1
set mode=$2
set period=$3

rm logs/hadd_check_bkg_${input}_mode_${mode}_$period
hadd -f output/${input}.${period}.mode_${mode}_check_bkg.root output_batch/${input}_${mode}/${period}/file/${input}.${period}.mode_${mode}_check_bkg.* > & logs/hadd_check_bkg_${input}_mode_${mode}_$period

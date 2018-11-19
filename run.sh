#/bin/csh

set input=$1
set mode=$2
set period=$3

rm logs/log_${input}_mode_${mode}_$period
./build_osc_ntuple $period $input $mode Card/18a_$period.card 0 1 0 >&logs/log_${input}_mode_${mode}_$period 

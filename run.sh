#/bin/csh

set input=$1
set mode=$2

rm logs/log_${input}_mode_$mode
./build_osc_ntuple sk4 $input $mode 18a_sk4.card 0 1 0 >&logs/log_${input}_mode_$mode 

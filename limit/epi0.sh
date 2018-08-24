#! /bin/bash

# The executable
LIMITWIER=./calclimit

# library
#LD_LIBRARY_PATH=/home/jsjang/gsl/lib:$LD_LIBRARY_PATH
LD_LIBRARY_PATH=/net/sukatmd1/work21/nishino/gsl/lib:$LD_LIBRARY_PATH

# Number of candidate events
nevents=0

# Rate, probably same for all limits
#min_rate=0.0
#max_rate=4.0
   
# Efficiency
#min_eff=0.0
#max_eff=1.0
#mean_eff=0.44
#sys_eff=0.18
#
# combine SK-I and SK-II
#(45.0*1489.2+43.6*798.6)/(1489.2+798.6)=44.51%
mean_eff=0.4451
sys_eff=0.19

# Exposure, probably same for all limits
#min_exp=7.0
#max_exp=10.0
#mean_exp=10.96
#sys_exp=0.01
#
# combine SK-I and SK-II
# (1489.2+798.6)/365.25*22.45*10^9*10/18*6.02*10^23 = 47.03
mean_exp=94.
#mean_exp=47.03
sys_exp=0.01

# Background
#min_bkg=0.0
#max_bkg=2.0
#num_bkg=465.626
num_bkg=4
oversamp_bkg=40
sys_bkg=1

$LIMITWIER \
        $nevents \
        $mean_eff $sys_eff \
        $mean_exp $sys_exp \
        $num_bkg $oversamp_bkg $sys_bkg 


# result
#
# max_rate = 2
# 9.115x10^33 yrs (no sys error)
# 8.380x10^33 yrs (w/ sys error)

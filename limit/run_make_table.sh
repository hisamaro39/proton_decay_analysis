#! /bin/bash

# The executable
LIMITWIER=./make_table

# Number of candidate events
nevents=0
   
# Efficiency
mean_eff=0.2
#systematic eff 
sys_eff=0.1

# Exposure, probably same for all limits
#sk1:91.7 sk2:49.2 sk3:31.9 sk4:133.5
mean_exp=133.5
#systematic eff
sys_exp=0.01

# Background
#unnormalized number of background
num_bkg=4 #default 4
#oversampling factor: sk4:63.7
oversamp_bkg=63.7
#systematic eff
sys_bkg=0.

# Use Baisian or Poisson (this is not working)
# Please comment-out or not POISSONLIMIT int code directory 
use_baisian=1

$LIMITWIER \
        $nevents \
        $mean_eff $sys_eff \
        $mean_exp $sys_exp \
        $num_bkg $oversamp_bkg $sys_bkg \
        $use_baisian


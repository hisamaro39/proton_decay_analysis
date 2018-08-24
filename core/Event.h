#ifndef _EVENT_
#define _EVENT_ 1

#include "Event.h"

/* 
	This Event.h is Specifically for the use with files created by 
	fillnt from the osc3d/ntuple in the repository. It represents
	a compressed set of data/mc used in the 3Flavor05 analysis
*/

struct Event{

/* Data Block */
           int		ipnu;
           float        dirnu[3];
           float        pnu;
           int          mode;
           int          ip;
           float        dprob;
           float        dir[3];
           float        amom;
           float        path;
           float        wall;
           int          itype;
           float        flxg[3];
           float        flxgo[3];
           float        flxh[3];
           float        flxho[3];
           float        weightx;
};


struct Nuclear{
    
          float   QSquared;
          int     NuclearEffect;
          float   Evis;
          float   NuclearEnergy;
}; 




#endif

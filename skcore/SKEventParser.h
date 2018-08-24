#ifndef _SKEventParser_
#define _SKEventParser_

#include <math.h>
#include <iostream>

#include "core/Event.h"
#include "core/EventParser.h"
#include "skcore/global.h"
#include "TChain.h"

#include "tools/DataManager.h"
#include "tools/TypeWrapper.h"

class SKEventParser : public EventParser
{
   public:
      SKEventParser() { ParserType = global::skP ;}
     ~SKEventParser() {;}

      virtual void   LoadDataManager(DataManager * );
		
      double GetEnergy()    { return (double) Energy;};
      double GetLeptonP()   { return (double) LeptonP;};
      double GetCosineZ()   { return (double) CosineTheta;};
      double GetNuCosineZ() { return (double) NuCosineTheta;};
      double GetHondaFluxRatio( int );
      double GetNeutrinoDir( int i ) { return (double) NuCosineTheta;}
      double GetEAveRMS()   { return (double) ErmsHax(0) ; }
      int    GetNEAve()     { return         nEAveHax(0) ; }

      double GetNeighborE   (int i) { return (double) NeighborE(i)    ;  }
      double GetNeighborEHax(int i) { return (double) NeighborEHax(i) ;  }

      int GetType() 	    { return EventType; };

      virtual void Parse();
      virtual void GetEventBinValues( std::vector<double> & vars );

   protected:

      float Energy;
      float LeptonP;

      float CosineTheta;
      float NuCosineTheta;

      TypeWrapper<int>   ipnu;
      TypeWrapper<int>   mode;
      TypeWrapper<int>   ip;
      TypeWrapper<int>   itype;

      TypeWrapper<float> dirnu;
      TypeWrapper<float> pnu;
      TypeWrapper<float> dir;
      TypeWrapper<float> amom;
      TypeWrapper<float> flxho;
      TypeWrapper<float> weightx;


      // for energy averaging 
      // must be build from FriendTree

      // The haxxed variables have been computed 
      // after introducing cuts on the neighbor energy 
      // distribution to use only neighbors between [0.5, 2.0]xParentE
      TypeWrapper<float> ErmsHax;
      TypeWrapper<int>   nEAveHax;

      TypeWrapper<float> Erms ;
      TypeWrapper<int>   nEAve;

      TypeWrapper<float> NeighborEHax;
      TypeWrapper<float> NeighborE;

      TypeWrapper<float> fsiweight;
      bool kUseFSIWeight;
};


#endif

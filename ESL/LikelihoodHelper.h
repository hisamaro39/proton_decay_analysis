#ifndef __LikelihoodHelper__
#define __LikelihoodHelper__

#include "tools/DataManager.h"

// Simple interface class for common routines 
// and methods used during the construction 
// and use of likelihood functions. 
// eg counting decay electrons

class LikelihoodHelper
{
 public:
   LikelihoodHelper( DataManager * , int );

   void LoadWrappers();

   int  muedcyBuild(  );
   int  mipBuild();
   bool FCFVCut();
   bool SingleRingCut();
   bool SubGeVCut();
   
   //
   float GetEmax ()       { return emax; }
   float GetTotalEnergy() { return tot_energy;}
   int   GetMERIndex()    { return mering;}
 

 protected:

   DataManager * dm;

   int   skgen;
   float emax;
   float tot_energy;
   int   mering;

   TypeWrapper<Int_t>    nmue; 
   TypeWrapper<Int_t>    nmue_sel; 
   TypeWrapper<Int_t>    nring; 
   TypeWrapper<Float_t>  evis; 
   TypeWrapper<Float_t>  etime; 
   TypeWrapper<UInt_t>   etype;
   TypeWrapper<Float_t>  ehit;
   TypeWrapper<Float_t>  egood;

   TypeWrapper<Float_t> wall;
   TypeWrapper<UInt_t>  nhitac;
   TypeWrapper<UInt_t>   ip;
   TypeWrapper<Float_t>  amome;
   TypeWrapper<Float_t> prmslg; 

};


#endif

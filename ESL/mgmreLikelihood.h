#ifndef _mgmreLikelihood_
#define _mgmreLikelihood_

#include <iostream>
#include <math.h>

#include "ESL/LikelihoodReturn.h"
#include "ESL/LikelihoodHelper.h"

#include "tools/TypeWrapper.h"


class mgmreLikelihood : public LikelihoodReturn
{

  public:

    // for likelihood construction
    mgmreLikelihood(DataManager *, int );

    // initialize using the default likelihood variable in the cardreader 
    mgmreLikelihood(DataManager *, CardReader *, int );
 
    // to specify which variable name the likilehood file should
    // be read in from
    mgmreLikelihood(DataManager *, CardReader *, int , const char * );
   
    int   GetDecayE() { return muedcy; }
    int   GetMER() { return mering; }
       
    // interface with the DataManager
    void LoadWrappers(); 
    virtual float llBuild(  );
    
    //functions used to compute input variables to the likelihood
    int   mipBuild();
    float mpidBuild();
    float fmomBuild();
    float dposeBuild();
    float pteBuild();
    
    int ringnBuild();
    float transmomBuild();

    int   GetEnergyBin( );

    // variables with persistent and multiple use in the 
    // code body
    float emax;
    float tot_energy;
    int   mering;
    int   muedcy;
    int   skgen;


    float   llmpid   ; 
    float   llmue    ; 
    float   llfmom   ; 
    float   lldpos   ; 

    float   GetLLmpid() { return llmpid   ;  }
    float   GetLLmue () { return llmue    ;  }
    float   GetLLfmom() { return llfmom   ;  } 
    float   GetLLdpos() { return lldpos   ;  }

    void    Init();
    void    FillGlobals();

    // For construction of the likelihood histograms
    bool  Precuts();
    int   IsSignalBG();

    void  DefineHistograms() ;
    void  FillHistograms() ;

    // Wrappers to variables from data manager
    TypeWrapper<Int_t>    nring; 
    TypeWrapper<Float_t>  amome;
    TypeWrapper<Int_t>    nmue; 
    TypeWrapper<Int_t>    nmue_sel; 
    TypeWrapper<Float_t>  evis; 
    TypeWrapper<Float_t>  etime; 
    TypeWrapper<UInt_t>   etype;
    TypeWrapper<Float_t>  ehit;
    TypeWrapper<Float_t>  egood;
    TypeWrapper<Float_t>  pos;
    TypeWrapper<Float_t>  dirtotmue;
    TypeWrapper<Float_t>  etotmue;
    TypeWrapper<Float_t>  amomm; 

    TypeWrapper<Float_t> msdir;
    TypeWrapper<Float_t> prmslg; 
    TypeWrapper<Float_t> probms; 
    TypeWrapper<Float_t> epos; 

    TypeWrapper<Float_t> dir;
    
    TypeWrapper<Float_t> oscwgt;
    TypeWrapper<Int_t>   ipnu;
    TypeWrapper<Int_t>   mode;

    LikelihoodHelper * llh;

 protected:
    int nEnergyBins;
    int nHistograms;  //4
    int sigbgbase[2];

    int        baseList[4]  ;
    const char * Titles[4]  ;
    const char * varlist[4] ;
    const char * SB[2]      ;
};

#endif

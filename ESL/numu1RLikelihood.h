#ifndef _numu1RLikelihood_
#define _numu1RLikelihood_

#include <iostream>
#include <math.h>

#include "ESL/LikelihoodReturn.h"

#include "tools/TypeWrapper.h"
#include "TFile.h"


class numu1RLikelihood : public LikelihoodReturn
{

  public:

    // initialize using the default likelihood variable in the cardreader 
    numu1RLikelihood(DataManager *, CardReader *, int );
 
    // to specify which variable name the likilehood file should
    // be read in from
    numu1RLikelihood(DataManager *, CardReader *, int , const char * );
   
       
    // interface with the DataManager
    void LoadWrappers(); 
    virtual float llBuild_numu1R();
    virtual void llBuild_numu1R_read();
    
    //functions used to compute input variables to the likelihood
    int   muedcyBuild();
    int   ndcye_muBuild();
    int   ndcye_piBuild();
    float etimeBuild();
    
    int   GetEnergyBin( );

    // variables with persistent and multiple use in the 
    // code body
    int   muedcy;
    int   skgen;

    // Wrappers to variables from data manager
    TypeWrapper<UInt_t>   ip;
    TypeWrapper<Int_t>    nring; 
    TypeWrapper<Float_t>  amome;
    TypeWrapper<Int_t>    nmue; 
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
    
    TypeWrapper<Int_t>   ipnu;
    TypeWrapper<Int_t>   mode;
    
    TypeWrapper<Float_t> oscwgt;

    
};

#endif

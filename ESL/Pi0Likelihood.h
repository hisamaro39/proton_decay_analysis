#ifndef _Pi0Likelihood_
#define _Pi0Likelihood_

#include <iostream>
#include <math.h>

#include "ESL/LikelihoodReturn.h"
#include "ESL/LikelihoodHelper.h"
#include "tools/TypeWrapper.h"


class Pi0Likelihood : public LikelihoodReturn
{

  public:
    
    // store skgeneration for access to LikelihoodHelper
    Pi0Likelihood(DataManager * , int );

    Pi0Likelihood(DataManager * , CardReader * );

    // indicate which variable contains the likelihood filename 
    Pi0Likelihood(DataManager *, CardReader * , const char *);

    int Pi0Selection();

    // interface with the DataManager
    void LoadWrappers(); 
    virtual float llBuild();

    //functions used to compute input variables to the likelihood
    float ring2fracBuild();
    float pi0massBuild();
    float deltapolfitBuild();

    int   GetEnergyBin( );

    bool  Precuts();
    int   IsSignalBG();

    void  DefineHistograms() ;
    void  FillHistograms() ;
 
 protected:

    // Wrappers to variables from data manager
    TypeWrapper<Float_t>  amome;
    TypeWrapper<Float_t>  pi0like;
    TypeWrapper<Float_t>  pi0_e;
    TypeWrapper<Float_t>  pi0mass;

    TypeWrapper<Int_t>   ipnu;
    TypeWrapper<Int_t>   mode;
    TypeWrapper<Float_t> evis;
    TypeWrapper<Float_t> wall;
    TypeWrapper<Int_t>   nring;
    TypeWrapper<UInt_t>  nhitac;
    TypeWrapper<UInt_t>  ip; 

    LikelihoodHelper * llh;
};

#endif

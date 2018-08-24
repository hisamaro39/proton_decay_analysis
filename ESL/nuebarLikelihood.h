#ifndef _nuebarLikelihood_
#define _nuebarLikelihood_

#include <iostream>
#include <math.h>

#include "ESL/LikelihoodReturn.h"
#include "ESL/LikelihoodHelper.h"
#include "ESL/mgmreLikelihood.h"
#include "tools/TypeWrapper.h"


class nuebarLikelihood : public LikelihoodReturn
{

  public:
    // for constructing the likelihood file
    nuebarLikelihood( DataManager * _dm , int skgen, const char *  );

    //
    nuebarLikelihood(DataManager *, CardReader * , int);
   
    // indicate which variable contains the likelihood filename 
    nuebarLikelihood(DataManager *, CardReader * , int,  const char *);


    // interface with the DataManager
    void LoadWrappers(); 
    virtual float llBuild_nue();
    virtual float llBuild_nuebar();

    void    Init();
    void    FillGlobals();

    // For construction of the likelihood histograms
    bool  Precuts();
    int   IsSignalBG();

    void  DefineHistograms() ;
    void  FillHistograms() ;


    //functions used to compute input variables to the likelihood
    int   ringnBuild();
    float transmomBuild();

    int   GetEnergyBin( );
    
    // variables with persistent and multiple use in the 
    // code body
    float emax;
    float tot_energy;
    int   mering;
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
    TypeWrapper<Float_t> oscwgt;

    TypeWrapper<Int_t>   ipnu;
    TypeWrapper<Int_t>   mode;

    LikelihoodHelper * llh;
    mgmreLikelihood  * llmre;
 protected:
    int nEnergyBins;
    int sigbgbase[4];
    const char * SB[4]       ;

    int nHistograms;  //3
    int          baseList[3] ;
    const char * Titles[3]   ;
    const char * varlist[3]  ;

};

#endif

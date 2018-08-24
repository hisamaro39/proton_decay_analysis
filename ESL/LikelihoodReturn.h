#ifndef _LikelihoodReturn_
#define _LikelihoodReturn_

#include <iostream>
#include <math.h>

#include "tools/CardReader.h"
#include "tools/DataManager.h"

#include "TFile.h"
#include "TH1.h"
#include "ESL/llHistogram.h"

/////////////////////////////
//
//  Parent class for returning the value of 
//  a likelihood function based on histograms
//  stored in a ROOT file
//
////////////////////////////////////



class LikelihoodReturn
{

  public:
    LikelihoodReturn(DataManager *, CardReader *);
   
    virtual float llBuild( );
    
    void  SetLikelihoodFile( TFile * a );
    void  SetOutputFile    ( const char * );
    void  SetBgOffset( int a ){ bgOffset = a; }
    void  SetSigOffset( int a ){ sigOffset = a; }
    void  SetVerbosity( bool x ){ kVerbose = x; }

    virtual void  ConstructLikelihood();
    virtual void  WriteHistograms();
    virtual void  DefineHistograms() {;}
    virtual void  FillHistograms() {;}
    virtual int   IsSignalBG() {return 0;}
    virtual bool  Precuts() {return 0;}

  protected:
    CardReader * read;
    DataManager * dm;

    TFile * llfile;
    TFile * outfile;
    std::map< int , TH1D* > h;
    std::vector< llHistogram* > vllH;

    float LoadLikelihood(int Prefix, float SearchKey, int eBin);    
    
    bool kVerbose;  
    int  bgOffset;
    int  sigOffset;

};

#endif

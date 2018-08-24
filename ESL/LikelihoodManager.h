#ifndef _LikelihoodManager_
#define _LikelihoodManager_

#include "ESL/mgmreLikelihood.h"
#include "ESL/Pi0Likelihood.h"
#include "ESL/nuebarLikelihood.h"
#include "ESL/numuLikelihood.h"
#include "ESL/numu1RLikelihood.h"
#include "ESL/pcstopLikelihood.h"
#include "ESL/pcthruLikelihood.h"

#include "tools/DataManager.h"
#include "tools/CardReader.h"
#include <sstream>
#include <string>


class LikelihoodManager
{
  public:

    // int is the SK generation [1,4]
    LikelihoodManager( int );

    // int is the SK generation [1,4]
    LikelihoodManager( const char * , int );
    LikelihoodManager( CardReader * , int );
   ~LikelihoodManager(     );

    void SetDataManager( DataManager * _dm ) { dm = _dm; Initialize(); }


    // MultiRing MultiGeV information
    float GetMGMRELikelihood();
    float GetMGMRELikelihood_nue();
    float GetMGMRELikelihood_nuebar();
    float GetnumuLikelihood_numu(); 
    float Getnumu1RLikelihood();
    float GetpcstopLikelihood();
    float GetpcthruLikelihood();
    int   GetDecayE(){ return nDecayE; }

    //
    // Pi0 Like  information
    //
   
    // Whether or not this event is selected as a Pi0-like
    int   GetPi0Selection();

    // The value of the Pi0 likelihood 
    float   GetPi0Likelihood();

    float   llmpid   ; 
    float   llmue    ; 
    float   llfmom   ; 
    float   lldpos   ; 

    float   GetLLmpid() { return llmpid   ;  }
    float   GetLLmue () { return llmue    ;  }
    float   GetLLfmom() { return llfmom   ;  } 
    float   GetLLdpos() { return lldpos   ;  }

  private:
    void Initialize();

    bool kUseMuonLikelihoods;

    CardReader      * cr;

    DataManager     * dm;

    mgmreLikelihood  * llmgmre; 
    Pi0Likelihood    * llpi0;
    nuebarLikelihood * llnue;

    // These are curently unused in oscillation
    // analysis 
    numuLikelihood   * llnumu;
    numu1RLikelihood * llnumu1R;
    pcstopLikelihood * llpcstop;
    pcthruLikelihood * llpcthru;

    int skx;
    int nDecayE;

    bool kInitialized;

    std::stringstream ss; 
};
#endif


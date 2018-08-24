#include "ESL/LikelihoodManager.h"
#include <stdlib.h>

LikelihoodManager::LikelihoodManager( int skg )
{

   ss << getenv("ATMPD_ROOT");
   ss << "/src/analysis/EventSelectionLikelihoods/Cards/Default.card\0" ;

   cr = new CardReader( ss.str().c_str() );
   skx = skg -1 ;
   llmgmre = 0;
   llpi0 = 0;
   
   llnue = 0;
   llnumu = 0;
   llnumu1R = 0; 
   llpcstop = 0;
   llpcthru = 0;
   
   kInitialized = false;
}


LikelihoodManager::LikelihoodManager( const char * filename , int skg )
{
   cr = new CardReader( filename );
   skx = skg -1 ;
   llmgmre = 0;
   llpi0 = 0;
   llnue = 0;
   llnumu=0;
   llnumu1R = 0; 
   llpcstop = 0;   
   llpcthru = 0;
   
   kInitialized = false;
}


LikelihoodManager::LikelihoodManager( CardReader * _cr , int skg )
{
   cr = _cr;
   skx = skg -1 ;
   llmgmre = 0;
   llpi0 = 0;
   llnue = 0;
   llnumu=0;
   llnumu1R = 0; 
   llpcstop = 0;   
   llpcthru = 0;
   
   kInitialized = false;
}


LikelihoodManager::~LikelihoodManager()
{
   delete cr;

   if( llmgmre ) delete llmgmre; 
   if( llpi0   ) delete llpi0; 
   if( llnue   ) delete llnue;
   if( llnumu  ) delete llnumu;
   if( llnumu1R) delete llnumu1R;
   if( llpcstop) delete llpcstop;
   if( llpcthru) delete llpcthru;
}


void LikelihoodManager::Initialize()
{
   if( llmgmre ) delete llmgmre; 
   if( llpi0   ) delete llpi0; 
   if( llnue   ) delete llnue;
   if( llnumu  ) delete llnumu;
   if( llnumu1R) delete llnumu1R;
   if( llpcstop) delete llpcstop; 
   if( llpcthru) delete llpcthru;  

   llpi0    = new Pi0Likelihood   ( dm , cr ,       "Pi0Likelihood"   );
   llmgmre  = new mgmreLikelihood ( dm , cr , skx , "mgmreLikelihood" ); 
   llnue    = new nuebarLikelihood( dm , cr , skx , "nuebarLikelihood");

   kUseMuonLikelihoods = false;
   cr->GetKey("UseMuonLikelihoods" ,  kUseMuonLikelihoods ); 

   if( kUseMuonLikelihoods )
   {
     llnumu   = new numuLikelihood  ( dm , cr , skx , "numuLikelihood"  );
     llnumu1R = new numu1RLikelihood( dm , cr , skx , "numu1RLikelihood");
     llpcstop = new pcstopLikelihood( dm , cr , skx , "pcstopLikelihood");
     llpcthru = new pcthruLikelihood( dm , cr , skx , "pcthruLikelihood");
   }

   kInitialized = true;
}


float LikelihoodManager::GetMGMRELikelihood()
{
    if( !kInitialized )
    {
       std::cout << "Warning! LikelihoodManager::GetMGMRELikelihood , DataManager"
                 << " has not been initialized and likelihood histograms are unloaded."
                 << " returning 0.  Please be sure to call LikelihoodManager::SetDataManager() "
                 << std::endl;
       return 0.;
    }


    float ll = llmgmre->llBuild();
    nDecayE  = llmgmre->GetDecayE();

    llmpid  = llmgmre->GetLLmpid() ; 
    llmue   = llmgmre->GetLLmue () ; 
    llfmom  = llmgmre->GetLLfmom() ; 
    lldpos  = llmgmre->GetLLdpos() ; 



    return ll;
}

float LikelihoodManager::GetMGMRELikelihood_nue()
{
    if( !kInitialized )
    {
       std::cout << "Warning! LikelihoodManager::GetMGMRELikelihood_nue , DataManager"
                 << " has not been initialized and likelihood histograms are unloaded."
                 << " returning 0.  Please be sure to call LikelihoodManager::SetDataManager() "
                 << std::endl;
       return 0.;
    }


    float ll = llnue->llBuild_nue();

    return ll;
}

float LikelihoodManager::GetMGMRELikelihood_nuebar()
{
    if( !kInitialized )
    {
       std::cout << "Warning! LikelihoodManager::GetMGMRELikelihood_nue , DataManager"
                 << " has not been initialized and likelihood histograms are unloaded."
                 << " returning 0.  Please be sure to call LikelihoodManager::SetDataManager() "
                 << std::endl;
       return 0.;
    }


    float ll = llnue->llBuild_nuebar();
    return ll;
}


float LikelihoodManager::GetnumuLikelihood_numu()
{
    if( !kInitialized || !kUseMuonLikelihoods )
    {
       std::cout << "Warning! LikelihoodManager::GetnumuLikelihood_numu , DataManager"
                 << " has not been initialized and likelihood histograms are unloaded."
                 << " returning 0.  Please be sure to call LikelihoodManager::SetDataManager() "
                 << std::endl;
       return 0.;
    }
 
 
    float ll = llnumu->llBuild_numu();
    
    return ll;

}


float LikelihoodManager::Getnumu1RLikelihood()
{
    if( !kInitialized || !kUseMuonLikelihoods )
    {
       std::cout << "Warning! LikelihoodManager::GetnumuLikelihood_numu , DataManager"
                 << " has not been initialized and likelihood histograms are unloaded."
                 << " returning 0.  Please be sure to call LikelihoodManager::SetDataManager() "
                 << std::endl;
       return 0.;
    }
 
 
    float ll = llnumu1R->llBuild_numu1R();
    
    return ll;

}


float LikelihoodManager::GetpcstopLikelihood()
{
    if( !kInitialized || !kUseMuonLikelihoods )
    {
       std::cout << "Warning! LikelihoodManager::GetnumuLikelihood_numu , DataManager"
                 << " has not been initialized and likelihood histograms are unloaded."
                 << " returning 0.  Please be sure to call LikelihoodManager::SetDataManager() "
                 << std::endl;
       return 0.;
    }
    
    float ll = llpcstop->llBuild_pcstop();
    
    return ll;
                 
}



float LikelihoodManager::GetpcthruLikelihood()
{
    if( !kInitialized || !kUseMuonLikelihoods )
    {
       std::cout << "Warning! LikelihoodManager::GetnumuLikelihood_numu , DataManager"
                 << " has not been initialized and likelihood histograms are unloaded."
                 << " returning 0.  Please be sure to call LikelihoodManager::SetDataManager() "
                 << std::endl;
       return 0.;
    }
    
    float ll = llpcthru->llBuild_pcthru();
    
    return ll;
                 
}


int LikelihoodManager::GetPi0Selection()
{
    if( !kInitialized )
    {
       std::cout << "Warning! LikelihoodManager::GetPi0Selection , DataManager"
                 << " has not been initialized and likelihood histograms are unloaded."
                 << " returning 0.  Please be sure to call LikelihoodManager::SetDataManager() "
                 << std::endl;
       return 0;
    }
    return llpi0->Pi0Selection();
}

float LikelihoodManager::GetPi0Likelihood()
{
    if( !kInitialized )
    {
       std::cout << "Warning! LikelihoodManager::GetPi0Likelihood , DataManager"
                 << " has not been initialized and likelihood histograms are unloaded."
                 << " returning 0.  Please be sure to call LikelihoodManager::SetDataManager() "
                 << std::endl;
       return 0.;
    }
    return llpi0->llBuild();
}

#include "ESL/Pi0Likelihood.h"

Pi0Likelihood::Pi0Likelihood( DataManager * _dm , int skgen ):LikelihoodReturn( _dm, 0 )
{
   LoadWrappers();
   std::cout << "Pi0Likelihood:: Initializing for likelihood construction " <<std::endl;
   llh = new LikelihoodHelper( dm , skgen );
}


Pi0Likelihood::Pi0Likelihood( DataManager * _dm , CardReader * reader  )
:LikelihoodReturn( _dm, reader )
{
  
   LoadWrappers();
   std::string input;

   read->GetKey("Pi0Likelihood", input );
   std::cout << "Pi0Likelihood:: Using likelihood information from" <<std::endl;
   std::cout << "                  " << input << std::endl;

   TFile * file = new TFile( input.c_str()  ); 
   SetLikelihoodFile( file ); 
}


Pi0Likelihood::Pi0Likelihood( DataManager * _dm , CardReader * reader , const char * varname )
:LikelihoodReturn( _dm, reader )
{
   LoadWrappers();
  
   std::string input;

   read->GetKey(varname, input );
   std::cout << "Pi0Likelihood:: Using likelihood information from" <<std::endl;
   std::cout << "                  " << input << std::endl;

   TFile * file = new TFile( input.c_str()  ); 
   SetLikelihoodFile( file ); 
}


// requires a DataManager to be defined
void Pi0Likelihood::LoadWrappers()
{
   dm->Get("amome"  , amome   );
   dm->Get("pi0like", pi0like );
   dm->Get("pi0mass", pi0mass );
   dm->Get("pi0_e"  , pi0_e   );

   dm->Get("mode"  , mode   );
   dm->Get("ipnu"  , ipnu   );
   dm->Get("ip"  , ip   );

// kVerbose = true ;

}

int  Pi0Likelihood::Pi0Selection()
{

   int Selection = 0;

   // EnergyBin returns histogram bin numbers
   //  71-75, so chop off the 70
   int bin  = GetEnergyBin() - 70;
   float ll = llBuild(); 

   if( pi0mass(0) <= 100. ) return 0;
   
   // the following cases are orthogonal 
   // and their respective outputs are the same
   if( bin == 1               ) Selection = 1; 
   if( bin == 2 && ll <  0.   ) Selection = 1; 
   if( bin == 3 && ll < -0.5  ) Selection = 1; 
   if( bin >= 4 && ll < -1.0  ) Selection = 1; 

   if( kVerbose && Selection == 1 )
      std::cout << "Pi0Likelihood::Pi0Selection " 
                << log10(amome(0)) << " " << pi0mass(0) << " " << ll << " " << bin << std::endl;


   // otherwise we are not pi0;
   return Selection;
}




float Pi0Likelihood::llBuild( )
{
    float ll = 0.;

    float ring2frac   = ring2fracBuild();
    float pi0mass     = pi0massBuild();
    float deltapolfit = deltapolfitBuild();

    int bin = GetEnergyBin( );
    SetBgOffset( 100 );    

    // this is the actual likelihood value computation
    if( bin > 0 )
    {
       ll  = 0.;
       ll += LoadLikelihood( 6100 , ring2frac     , bin );
       ll += LoadLikelihood( 7100 , pi0mass       , bin );
       ll += LoadLikelihood( 8100 , deltapolfit   , bin );
    }

   return ll;
}


float Pi0Likelihood::ring2fracBuild()
{

   float fraction  = pi0_e(0,1);
   fraction /= ( pi0_e(0,0) + pi0_e(0,1) );


   return fraction;
}



float Pi0Likelihood::pi0massBuild()
{
   float pi0 = pi0mass(0);

   if( pi0 >= 300.0 ) 
      pi0 = 299.0; 

   return pi0;

}



float Pi0Likelihood::deltapolfitBuild()
{

   float delta = pi0like(0)-pi0like(1);

   if ( delta >=  400.0 ) delta =  399.0;
   if ( delta <= -200.0 ) delta = -199.0;

   return delta;
}




int Pi0Likelihood::GetEnergyBin( )
{

  float logE = log10( amome(0) );

  int nbins = 5;
  float edges[] ={ 0., 2.4, 2.6, 2.8, 3.0, 1.e6 };

  int bin = 0;

  for( int i = 1 ; i <= nbins ; i++ )
    if( logE > edges[i-1] && logE < edges[i] )
      bin = 70 + i;  // add 70 to match histogram 
                     // numbering

  // adjust to match histogram numbering
  return bin;

}

bool Pi0Likelihood::Precuts()
{
   bool kPass = false;

   if( llh->FCFVCut() )
    if( llh->SingleRingCut() )
     if( llh->SubGeVCut()   )
      if( llh->muedcyBuild()  == 0  )
       if( amome(0) > 100.0 && ip(0) == 2 )
        if( pi0mass(0) > 100.0 ) 
           kPass = true ;


   return kPass; 
}


void Pi0Likelihood::DefineHistograms()
{

   int nEnergyBins = 5;
   int ID;

   int nHistograms = 3;
   int baseList[]  = { 6000 , 7000 ,  8000 };
   const char * Titles  []  = {  "Second Ring Fraction"   , 
                                 "POLfit pi0 mass"        ,
                                 "Delta POLfit Likelihood" };

   const char * SB []  = { "Signal" , "Background" };
   int sigbgbase   []  = { 170      ,  270         };
   int nBinsList   []  = { 30    , 20    , 40     };
   float LowEdges  []  = { 0.0   , 100.0 , -200.0 };  
   float HiEdges   []  = { 0.7   , 200.0 ,  400.0 };  
                                   
   const char * varlist [] = { "ring2frac" , "pi0mass" , "deltapolfit" };

   llHistogram * ph;
   std::stringstream ss;
   for( int i = 1 ; i <= nEnergyBins ; i++ )
    for( int j = 0 ; j < nHistograms ; j ++ ) 
     for( unsigned k = 0 ; k < 2 ; k ++ ) 
     {
     
       ss.str(""); ss << Titles[j] << ", " << SB[k] << ", Bin " << i ;
       ID = baseList[j] + sigbgbase[k] + i ;
       ph = new llHistogram( ID, nBinsList[j], LowEdges[j], HiEdges[j], ss.str().c_str() , varlist[j] ); 
    
       ph->h->SetLineWidth(2);
       if ( k == 0 ) ph->h->SetLineColor( 2 );
       if ( k == 1 ) ph->h->SetLineColor( 4 );
      
       h[ID] = ph->h; 
       vllH.push_back( ph );
     }
  

}



void Pi0Likelihood::FillHistograms() 
{

   float ring2frac   = ring2fracBuild();
   float pi0mass     = pi0massBuild();
   float deltapolfit = deltapolfitBuild();
   int bin           = GetEnergyBin( );

   int baseList[]    = { 6000  , 7000 ,  8000 };
   int sigbgbase []  = { 100       , 200      };
   float lvarlist []  = { ring2frac , pi0mass , deltapolfit };

   int k = IsSignalBG();
   if( k < 0    ) return;
   if( bin <= 0 ) return;

   for ( unsigned j = 0 ; j < 3 ; j ++ )
   {
      int ID = baseList[j] + sigbgbase[k] + bin ;
      h [ ID ]->Fill( lvarlist[j] );
   }
  

}

//  signal : 0
//  BG : 1
//  other: -1
int Pi0Likelihood::IsSignalBG()
{
   // signal
   if( abs(ipnu(0)) == 12  && abs( mode(0) )  == 1 )  return 0 ;
  
   // BG
   if( abs(mode(0)) > 30 )  return 1 ;

   return -1; 
}







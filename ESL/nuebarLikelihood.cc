#include "ESL/nuebarLikelihood.h"
#include <fstream>

using namespace std;

nuebarLikelihood::nuebarLikelihood( DataManager * _dm , int skgen, const char * input ):LikelihoodReturn( _dm, 0 )
{
   LoadWrappers();
   std::cout << "nuebarLikelihood:: Initializing for likelihood construction " <<std::endl;
   llh = new LikelihoodHelper( dm , skgen );
   llmre = new mgmreLikelihood( dm , skgen );

   TFile * file = new TFile( input  ); 
   llmre->SetLikelihoodFile( file ); 
 
   Init();
}

nuebarLikelihood::nuebarLikelihood( DataManager * _dm , CardReader * reader , int skg )
:LikelihoodReturn( _dm, reader )
{
   skgen = skg;
   llh = new LikelihoodHelper( dm , skgen );
 
   LoadWrappers();
  
   std::string input;

   read->GetKey("nuebarLikelihood", input );
   std::cout << "nuebarLikelihood:: Using likelihood information from" <<std::endl;
   std::cout << "                  " << input << std::endl;

   TFile * file = new TFile( input.c_str()  ); 
   SetLikelihoodFile( file ); 
}


nuebarLikelihood::nuebarLikelihood( DataManager * _dm , CardReader * reader , int skg, const char * varname )
:LikelihoodReturn( _dm, reader )
{
   skgen = skg;
   llh = new LikelihoodHelper( dm , skgen );
   
   LoadWrappers();
  
   std::string input;

   read->GetKey(varname, input );
   read->GetKey("nuebarLikelihood", input );
   std::cout << "nuebarLikelihood:: Using likelihood information from" <<std::endl;
   std::cout << "                  " << input << std::endl;

   TFile * file = new TFile( input.c_str()  ); 
   SetLikelihoodFile( file ); 
}


void nuebarLikelihood::Init()
{

   nEnergyBins = 5;
   nHistograms = 3;

   baseList[0] = 180000;
   baseList[1] = 190000;
   baseList[2] = 200000;

   Titles[0] =  "Number of Decay E"   ; 
   Titles[1] =  "Number of Rings" ;
   Titles[2] =  "Transverse Momentum";

   varlist[0] = "muedcy"    ;
   varlist[1] = "nring"     ;
   varlist[2] = "transmom"  ;

   SB[0] = "CC Nue Signal"      ;
   SB[1] = "CC Nuebar Signal"  ;
   SB[2] = "Background"  ;
   SB[3] = "Unused"      ;

   sigbgbase[0]  =  0;
   sigbgbase[1]  = 10;
   sigbgbase[2]  = 20;
   sigbgbase[3]  = 50;

}


// requires a DataManager to be defined
void nuebarLikelihood::LoadWrappers()
{
   dm->Get("nmue", nmue);
   dm->Get("evis", evis);
   dm->Get("ip", ip);
   dm->Get("nring", nring);
   dm->Get("amome", amome);
   dm->Get("amomm", amomm);
   dm->Get("nmue", nmue);
   dm->Get("evis", evis);
   dm->Get("etime", etime);
   dm->Get("etype", etype);
   dm->Get("ehit", ehit);
   dm->Get("egood", egood);
   dm->Get("pos", pos);
   dm->Get("dirtotmue", dirtotmue);
   dm->Get("etotmue", etotmue);

   dm->Get("msdir"  , msdir   );
   dm->Get("prmslg" , prmslg  ); 
   dm->Get("epos"   , epos    );

   dm->Get("dir"    , dir     );
   dm->Get("oscwgt" , oscwgt  );

   dm->Get("mode"  , mode   );
   dm->Get("ipnu"  , ipnu   );

}



float nuebarLikelihood::llBuild_nue( )
{
    float ll = 0.;
    
   

    int   mip    = llh->mipBuild();
    FillGlobals();
          muedcy = llh->muedcyBuild();
    float ringn  = ringnBuild();
    float transmom = transmomBuild();

    sigOffset =  0 ; // cc nue 
    bgOffset  = 20 ; // BG

    int bin = GetEnergyBin( );
    
    // this is the actual likelihood value computation
    if( bin > 0 )
    {
       ll += LoadLikelihood( 180000 , muedcy     , bin );
       ll += LoadLikelihood( 190000 , ringn      , bin );
       ll += LoadLikelihood( 200000 , transmom   , bin );
    }

   return ll;
}

float nuebarLikelihood::llBuild_nuebar( )
{
    float ll = 0.;
   

    int   mip      = llh->mipBuild();
    FillGlobals();
          muedcy   = llh->muedcyBuild();
    float ringn    = ringnBuild();
    float transmom = transmomBuild();

    int bin = GetEnergyBin( ) ;

    sigOffset = 10 ; // nuebar CC 
    bgOffset  = 20 ; // BG
 
    // this is the actual likelihood value computation
    if( bin > 0 )
    {
       ll += LoadLikelihood( 180000 , muedcy     , bin );
       ll += LoadLikelihood( 190000 , ringn      , bin );
       ll += LoadLikelihood( 200000 , transmom   , bin );
    }

   return ll;
}




int nuebarLikelihood::GetEnergyBin( )
{

  int nbins = 5;
  float edges[] ={ 1.00, 2.5118, 5.0118, 10.00, 19.952, 1.0e6 };
  float E = etotmue(0)/1000.0;

  int bin = 0;

  for( int i = 1 ; i <= nbins ; i++ )
    if( E > edges[i-1] && E < edges[i] )
      bin = i;

  // adjust to match histogram numbering
  return bin;

}

//======================================================================================


int nuebarLikelihood::ringnBuild()
{
    int lringn;
    
    lringn = nring(0);
    
    return lringn;
}


float nuebarLikelihood::transmomBuild()
{
    int i;
    
    float cos_angle[10], sin_angle[10];
    float ltransmom=0;

    for(int j=0; j< nring(0); j++)
    {
        cos_angle[j]= dir(mering,0)*dir(j,0) + dir(mering,1)*dir(j,1) + dir(mering,2)*dir(j,2);
	
	if(cos_angle[j]>=1)
	{
	    cos_angle[j]=1;
	}
	
	sin_angle[j]= sqrt(1-cos_angle[j]*cos_angle[j]);
	
	ltransmom = ltransmom + amome(j)*sin_angle[j];
    } 
        
    return ltransmom/tot_energy;     
}

bool nuebarLikelihood::Precuts()
{
   bool kPass = false;

   //std::string elike_mr ("mip == 2 && wall > 200. && nring > 1 && evis > 1330. ");
   if ( llmre->Precuts() )
    if ( llmre->llBuild() > 0.0 ) 
       kPass = true ;

   return kPass; 
}


void nuebarLikelihood::DefineHistograms()
{

   int ID;

   // in order  of base :   18 ,    19 ,   20 ,  x 10000 
   int    nBinsList[]  = {  11  ,    7,     21 }; 
   double LowEdges []  = { -0.5 , -0.5, -0.025 }; 
   double HiEdges  []  = { 10.5 ,  6.5,  1.025 }; 

   llHistogram * ph;
   std::stringstream ss;
   
   for( unsigned i = 1 ; i <= nEnergyBins ; i++ )
    for( unsigned j = 0 ; j < nHistograms ; j ++ ) 
     for( unsigned k = 0 ; k < 4 ; k ++ ) 
     {
     
       ss.str(""); ss << Titles[j] << ", " << SB[k] << ", Bin " << i ;
       ID = baseList[j] + sigbgbase[k] + i ;
       ph = new llHistogram( ID, nBinsList[j], LowEdges[j], HiEdges[j], ss.str().c_str() , varlist[j] ); 
    
       ph->h->SetLineWidth(2);
       if ( k == 0 ) ph->h->SetLineColor( 2 );
       if ( k == 1 ) ph->h->SetLineColor( 4 );
       if ( k == 2 ) ph->h->SetLineColor( 3 );
      
       h[ID] = ph->h; 
       vllH.push_back( ph );
     }
  
}



void nuebarLikelihood::FillHistograms() 
{

   int   mip      = llh->mipBuild();
   FillGlobals();
         muedcy   = llh->muedcyBuild();
   float ringn    = ringnBuild();
   float transmom = transmomBuild();

   int   bin    = GetEnergyBin( );
   float lvars []  = { muedcy , ringn , transmom };

   int sigbgtype = IsSignalBG();
   if( sigbgtype < -5 ) return;
   if( bin <=  0      ) return;

   int ID;
   for ( unsigned j = 0 ; j < nHistograms ; j ++ )
   {
      //  nuesample
      if( sigbgtype == 1 )
         ID = baseList[j] + sigbgbase[0] + bin ;

      //  nuebarsample
      if( sigbgtype == -1 )
         ID = baseList[j] + sigbgbase[1] + bin ;

      // background CC numu || NC 
      if( abs ( sigbgtype) == 2 || sigbgtype == 0 )
         ID = baseList[j] + sigbgbase[2] + bin ;

      h [ ID ]->Fill( lvars[j] );
   }
  

}


int nuebarLikelihood::IsSignalBG()
{
   // signal
   if( ipnu(0) ==  12  && abs( mode(0) )  < 30 )   return  1 ;  // nue
   if( ipnu(0) == -12  && abs( mode(0) )  < 30 )   return -1 ;  // nuebar

   // background
   if( abs( mode(0) ) < 30 && ipnu(0) !=  12  ) return  2 ;  // bg to nue 
   if( abs( mode(0) ) < 30 && ipnu(0) != -12  ) return -2 ;  // bg to nuebar 
   if( abs( mode(0) ) > 30                    ) return  0 ;  // NC BG
 
   return -10;

}


void nuebarLikelihood::FillGlobals()
{
   tot_energy   = llh->GetTotalEnergy();
   mering       = llh->GetMERIndex();
}





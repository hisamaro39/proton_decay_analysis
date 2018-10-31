#include "ESL/mgmreLikelihood.h"

#include <fstream>

using namespace std;

mgmreLikelihood::mgmreLikelihood( DataManager * _dm , int skg ):LikelihoodReturn( _dm, 0 )
{
   skgen = skg;
   LoadWrappers();
   std::cout << "mgmreLikelihood:: Initializing for likelihood construction " <<std::endl;
   llh = new LikelihoodHelper( dm , skgen );

   Init();
}

mgmreLikelihood::mgmreLikelihood( DataManager * _dm , CardReader * reader , int skg )
:LikelihoodReturn( _dm, reader )
{
   skgen = skg;

   llh = new LikelihoodHelper( dm , skgen );

   LoadWrappers();
  
   std::string input;

   read->GetKey("mgmreLikelihood", input );
   std::cout << "mgmreLikelihood:: Using likelihood information from" <<std::endl;
   std::cout << "                  " << input << std::endl;

   TFile * file = new TFile( input.c_str()  ); 
   SetLikelihoodFile( file ); 
}



mgmreLikelihood::mgmreLikelihood( DataManager * _dm , CardReader * reader , int skg, const char * varname )
:LikelihoodReturn( _dm, reader )
{
   skgen = skg;
   llh = new LikelihoodHelper( dm , skgen );

   LoadWrappers();
  
   std::string input;

   read->GetKey( varname , input );
   std::cout << "mgmreLikelihood:: Using likelihood information from" <<std::endl;
   std::cout << "                  " << input << std::endl;

   TFile * file = new TFile( input.c_str()  ); 
   SetLikelihoodFile( file ); 
}

void mgmreLikelihood::Init()
{

   nEnergyBins = 5;
   nHistograms = 4;

   baseList[0] = 130000;
   baseList[1] = 140000;
   baseList[2] = 150000;
   baseList[3] = 160000;

   Titles[0] =  "MER PID"             ;
   Titles[1] =  "Number of Decay E"   ; 
   Titles[2] =  "MER Ring Fraction"   ; 
   Titles[3] =  "Max Decay E distance";

   varlist[0] = "prmslg" ;
   varlist[1] = "muedcy" ;
   varlist[2] = "fmom"   ;
   varlist[3] = "dpose"  ;

   SB[0] = "Signal"      ;
   SB[1] = "Background"  ;

   sigbgbase[0]  =  0;
   sigbgbase[1]  = 10;

}



// requires a DataManager to be defined
void mgmreLikelihood::LoadWrappers()
{
   dm->Get("nmue", nmue);
   dm->Get("evis", evis);
   dm->Get("nring", nring);
   dm->Get("amome", amome);
   dm->Get("amomm", amomm);
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



float mgmreLikelihood::llBuild( )
{
  std::cout << "mgmreLikelihood::llBuild" << std::endl;
    float ll = 0.;


//  std::cout << " --------------- crash: " << " mip " <<  std::endl;
    int   mip    = llh->mipBuild(); // must be called first!
    FillGlobals();
          
          muedcy = llh->muedcyBuild();  // muedcy is a class variable, and may be accessed externally
    float mpid   = mpidBuild();
    float dpose  = dposeBuild();
    float pte    = pteBuild();
    float fmom   = fmomBuild();
  
    sigOffset =  0 ; 
    bgOffset  = 10 ;
    
    int bin = GetEnergyBin( );

    // this is the actual likelihood value computation
    if( bin > 0 )
    {
       llmpid   = LoadLikelihood( 130000, mpid          , bin );
       llmue    = LoadLikelihood( 140000, float(muedcy) , bin );
       llfmom   = LoadLikelihood( 150000, fmom          , bin );

       lldpos   = 0.0;
       if( dpose > 1.0e-5 )
          lldpos  = LoadLikelihood( 150000, sqrt( dpose )  , bin );
 
       ll += LoadLikelihood( 130000, mpid          , bin );
       ll += LoadLikelihood( 140000, float(muedcy) , bin );
       ll += LoadLikelihood( 150000, fmom          , bin );
       
       if( dpose > 1.0e-5 )
          ll += LoadLikelihood( 160000, sqrt(dpose)   , bin );

      if( kVerbose )
         std::cout << "mgmre " << ll
                   << " " <<  LoadLikelihood( 130000, mpid          , bin ) 
                   << " " <<  LoadLikelihood( 140000, float(muedcy) , bin )
                   << " " <<  LoadLikelihood( 150000, fmom          , bin )
                   << " " <<  LoadLikelihood( 160000, sqrt(dpose)   , bin )
                   << std::endl;

    }

   return ll;
}


int mgmreLikelihood::GetEnergyBin( )
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


float mgmreLikelihood::mpidBuild()
{
   float fmom = sqrt(1.*prmslg(mering,1) ) - sqrt(1.*prmslg(mering,2) );

   return fmom;
}


float mgmreLikelihood::fmomBuild()
{
   // fractional momentum carried by mering
   float fmom = amome( mering ) / tot_energy;

   return fmom;
}


float mgmreLikelihood::dposeBuild()
{

//  return longest distance between primary vertex position
//  and decay electron position, divided by mering's energy
     
   float dpos;
   float dx;
   float x, y , z;
   int   nmuemax;

   int ehit_cut;


   if( skgen == 0 ) ehit_cut = 60;
   if( skgen == 1 ) ehit_cut = 30;
   if( skgen == 2 ) ehit_cut = 60;
   if( skgen == 3 ) ehit_cut = 60;


   nmuemax = ( nmue(0) > 10 ? 10 :  nmue(0) );
   if ( nmue(0) <= 0 ) return 0.;

   dpos=0.;
   for( int i = 0 ; i < nmuemax ; i++ )
   {
     if( skgen == 3) //for SK4 only
     {
         if( etime(i) < 0.1) continue;
	 if( etime(i) < 0.6) continue;
	 
	 if( etype(i) == 1 )
	 {
	    x = pos(0)-epos(i,0);
	    y = pos(1)-epos(i,1);
	    z = pos(2)-epos(i,2);
	    dx = sqrt(x*x + y*y + z*z);
	    
	    if( dx > dpos)
	    {
	        dpos = dx;
	    }
	 }
     }
     else
     {
   
     	 if( etime(i) <  0.1 ) continue;
     	 if( etime(i) >  0.8 && etime(i) < 1.2 ) continue;
     	 if( etype(i) == 1   && ehit(i) >= ehit_cut && egood(i) > 0.5 )
     	 {
         	x  = pos(0) - epos(i,0);
         	y  = pos(1) - epos(i,1);
         	z  = pos(2) - epos(i,2);
         	dx = sqrt( x*x + y*y + z*z );

		if( dx > dpos ) dpos = dx;
     	 }
     }
   }// end of loop on decay-e's


   return dpos / llh->GetEmax() ;
}


// actually not used to make the likelihood as of 20090416
float mgmreLikelihood::pteBuild()
{

  float pte=0.;
  float prod ;
 
  for( int i = 0 ; i <  nring(0) ; i++ )
  {
     prod = 0.;
     if( prmslg(i,1) > prmslg(i,2))     
     {
      for( int j = 0 ; j < 3 ; j ++ )
       prod += dirtotmue(j)*msdir(i,j,1);
   
      pte += amome(i)*sqrt(1.-prod*prod);

     }
     else
     {
      for( int j = 0 ; j < 3 ; j ++ )
       prod += dirtotmue(j)*msdir(i,j,2);

      pte += amomm(i)*sqrt(1.-prod*prod);
     }// end of e/mu- like testing
   }// end of ring loop

  return pte / pow( tot_energy, 0.75 );

}


bool mgmreLikelihood::Precuts()
{
   bool kPass = false;

   //std::string elike_mr ("mip == 2 && wall > 200. && nring > 1 && evis > 1330. ");

   if(   llh->FCFVCut()       )  
    if( ! llh->SingleRingCut() )  // MultiRings
     if( ! llh->SubGeVCut()     ) // MultiGeV
      if(   llh->mipBuild() == 2 )// Most energetic ring is e-like 
          kPass = true ;

   return kPass; 
}


void mgmreLikelihood::DefineHistograms()
{

   int ID;

   // in order  of base :   13 , 14 ,    15 ,     16  x 10000 
   int    nBinsList[]  = { 82  ,    11,     21,     84 }; 
   double LowEdges []  = { 0.0 , -0.5, -0.025, -0.025 }; 
   double HiEdges  []  = { 0.0 , 10.5,  1.025,  1.025 }; 

   // Binning changes as a function of energy for the PID variable
   double pid_low [5] = { -12.25, -20.5, -31.5,  -51.5, -101.5 };  
   double pid_high[5] = {   0.25,   0.5,   0.5,    2.5,    2.5 };  
                                   

   llHistogram * ph;
   std::stringstream ss;
   
   // process for all but PID variable first 
   for( unsigned i = 1 ; i <= nEnergyBins ; i++ )
    for( unsigned j = 1 ; j < nHistograms ; j ++ ) 
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
  
   // now build for the PRMSLG variable 
   unsigned j = 0 ;
   int nBins  = 82 ;
   for( unsigned i = 1 ; i <= nEnergyBins ; i++ )
     for( unsigned k = 0 ; k < 2 ; k ++ ) 
     {
     
       ss.str(""); ss << Titles[j] << ", " << SB[k] << ", Bin " << i ;
       ID = baseList[j] + sigbgbase[k] + i ;
       ph = new llHistogram( ID, nBins , pid_low[i-1] , pid_high[i-1], ss.str().c_str() , varlist[j] ); 
    
       ph->h->SetLineWidth(2);
       if ( k == 0 ) ph->h->SetLineColor( 2 );
       if ( k == 1 ) ph->h->SetLineColor( 4 );
      
       h[ID] = ph->h; 
       vllH.push_back( ph );
     }
  

}



void mgmreLikelihood::FillHistograms() 
{

   int   mip    = llh->mipBuild(); // must be called first!
   FillGlobals();
         muedcy = llh->muedcyBuild();  // muedcy is a class variable, and may be accessed externally
   float mpid   = mpidBuild();
   float dpose  = dposeBuild();
   float fmom   = fmomBuild();
   int   bin    = GetEnergyBin( );

   float lvars []  = { mpid , muedcy , fmom, sqrt(dpose) };

   int k = IsSignalBG();
   if( k < 0    ) return;
   if( bin <= 0 ) return;

   for ( unsigned j = 0 ; j < nHistograms ; j ++ )
   {
      int ID = baseList[j] + sigbgbase[k] + bin ;
      h [ ID ]->Fill( lvars[j] );
   }
  

}

//  signal : 0
//  BG : 1
//  other: -1
int mgmreLikelihood::IsSignalBG()
{
   // signal
   if( abs(ipnu(0)) == 12  && abs( mode(0) )  < 30 )  return 0 ;
   // BG
   else                                               return 1 ;
  
   return -1; 
}


void mgmreLikelihood::FillGlobals()
{
   tot_energy   = llh->GetTotalEnergy();
   mering       = llh->GetMERIndex();
   //cout << "mgmreLikelihood " << tot_energy << " " << mering << std::endl;

}





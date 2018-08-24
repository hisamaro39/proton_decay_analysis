#include "ESL/LikelihoodReturn.h"
#include <fstream>

using namespace std;


LikelihoodReturn::LikelihoodReturn( DataManager * _dm , CardReader * reader  )
{
   read = reader;
   dm = _dm;
   bgOffset = 10;
   sigOffset = 0;
// kVerbose = true;
   kVerbose = false;
   outfile  = 0 ;
}



float LikelihoodReturn::llBuild( )
{
   std::cout << "Warning LikelihoodReturn::llBuild , this function "
             << "must be overridden in a child class. Returning 0." << std::endl;

   return 0.;
}

float LikelihoodReturn::LoadLikelihood(int Prefix, float SearchKey, int eBin)
{
   float sContent;   // signal
   float bContent;   // background
   float sIntegral;
   float bIntegral;
   int   Bin;          // assumes binning for signal and bg is identical
   int   Signal = eBin + sigOffset ;
   int   Background = eBin + bgOffset; 

   // fetch the signal and Background histograms
   TH1D * _hS = h[ Prefix + Signal ];
   TH1D * _hB = h[ Prefix + Background ];
   
   Bin = _hS->FindBin( SearchKey );
   sContent = _hS->GetBinContent( Bin );
   bContent = _hB->GetBinContent( Bin );

   sIntegral = _hS->Integral();
   bIntegral = _hB->Integral();

   // we want the larger of 1.0e-0 and the Content for s or b
   if ( sContent < 1.0e-7 ) sContent = 1.0 ;
   if ( bContent < 1.0e-7 ) bContent = 1.0 ;

   if( Bin > _hS->GetNbinsX()  || Bin == 0 )
   {
     sContent = 1.0;
     bContent = 1.0;
   }
   
   if( kVerbose )
      std::cout << " Prefix: " << Prefix << std::endl
                << "   key: " << SearchKey <<  std::endl  
                << "   sHist: " << Signal << " " << " bHist: " << Background << std::endl  
                << "   sContent: " << sContent << "  bContent " << bContent << " "  
                << log10( sContent / sIntegral ) - log10( bContent / bIntegral ) << std::endl
                << "   Bin:  "  << Bin <<  "  nbins: " << _hS->GetNbinsX() << std::endl;
               

   return log10( sContent / sIntegral ) - log10( bContent / bIntegral );

}


#include "TKey.h"
void LikelihoodReturn::SetLikelihoodFile( TFile * a )
{
  std::cout << "SetLikelihoodFile" << std::endl;

  llfile = a;  
  int id;

  TIter  NextKey = llfile->GetListOfKeys();
  TKey * lKey;
  TObject* histogram;
  TH1D * h1d;

  std::string handle;
  std::string::size_type loc;
  // reload any xisting histograms
  h.clear();
  
  if( kVerbose )
     std::cout << " Now loading histograms from likelihood file " << llfile->GetName() << std::endl;
  // load in all of the histograms in the file
  // they should all be six digit numbers
  while(   lKey = (TKey*) NextKey() )
  { 
    histogram = lKey->ReadObj();
    if( histogram->IsA()->InheritsFrom( TH1::Class() ) )
    {
       h1d = (TH1D*) histogram; 
       handle = h1d->GetName();
      
       loc = handle.find("h");
       //remove any "h" added to histograms defined in 
       // paw and converted to root with h2root
       if( loc != std::string::npos )
         handle.erase( loc, 1 ); 
       id = atoi( handle.c_str() ); 

       // load the map with pointer to this histogram
       h[id] = h1d; 
      
       if( kVerbose )
          std::cout << "Loading : " << h1d->GetName() << " " << id << std::endl;
    }
  }

   return;
}


void LikelihoodReturn::ConstructLikelihood()
{

  DefineHistograms();
  
  for( unsigned i = 0 ; i < dm->GetEntries() ; i++ )
  {
    dm->GetEntry(i);
    if( ! Precuts() ) continue ; 
    FillHistograms( ); 
  }

}

void LikelihoodReturn::SetOutputFile( const char * x )
{
   outfile = new TFile( x , "recreate" );
}


void LikelihoodReturn::WriteHistograms()
{
   
  if ( outfile == 0 )
  {
    std::cout << "LikelihoodReturn::WriteHistograms  - warning, output file has not been defined! " << std::endl;
    std::cout << "  Set with LikelihoodRetrun::SetOutputFile() " << std::endl;
    return;
  }

  std::map< int , TH1D* >::iterator _i;

  outfile->cd();
  for( _i = h.begin() ; _i != h.end() ; _i++ )
  {
    float integral = _i->second->Integral();
    _i->second->Scale( 1.0 / integral );
    _i->second->Write();
  }

  outfile->Close();

}






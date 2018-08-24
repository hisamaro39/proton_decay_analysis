#include "skcore/SKEventParser.h"
#include <math.h>

void SKEventParser::Parse()
{
   // Transfer the TypeWrapper information
   // a bit wasteful but allows us to simultaneously
   // check for events that should be skipped 
   // and perform other flagging operations 
   SkipFlag = false;
   PDG = 0;

   MCweight = weightx(0);
   Mode     = mode(0);

   // Parent Neutrino Parameters 
   NuCosineTheta = -1. * dirnu(2); // + downgoing , - upgoing 
   Energy        = pnu(0);
   PDG           = ipnu(0);

   // Lepton Parameters 
   LeptonP = amom(0);
   CosineTheta = -1.*dir(2);  // + downgoing , - upgoing 

   EventType = itype(0) - 1;

   if ( EventType < 0 || EventType > global::EndOfBinTypes )
      SkipFlag = true;

   // no NC tau's allowed!
   if ( abs(Mode) > 30  && abs(PDG) == 16 ) SkipFlag = true; 
 
   MonteCarlo  = true;
   if ( Energy < 1.0e-7 && PDG == 0 && Mode == 0 )
     MonteCarlo  = false;

}


double SKEventParser::GetHondaFluxRatio( int NuType )
{

   if ( NuType == 1 )   
      return (double) flxho(1) / flxho(0);

   if ( NuType == 2 || NuType == 3 )
      return (double) flxho(0) / flxho(1);

   std::cerr << "Returning 0 from SKEventParser::GetHondaFluxRatio for type " << NuType << std::endl;
   return 0;

}


void SKEventParser::LoadDataManager( DataManager * dm )
{

   dm->Get("ipnu"   , ipnu     );
   dm->Get("mode"   , mode     );
   dm->Get("ip"     , ip       );
   dm->Get("itype"  , itype    );

   dm->Get("dirnu"  , dirnu    );
   dm->Get("pnu"    , pnu      );
   dm->Get("dir"    , dir      );
   dm->Get("amom"   , amom     );
   dm->Get("flxho"  , flxho    );
   dm->Get("weightx", weightx  );


   // only if friends
   dm->Get("ErmsHax" , ErmsHax  );
   dm->Get("nEAveHax", nEAveHax );

   dm->Get("Erms"    , Erms  );
   dm->Get("nEAve"   , nEAve );

   dm->Get("NeighborEHax" , NeighborEHax );
   dm->Get("NeighborE"    , NeighborE    );

   Parse( );	
}


void SKEventParser::GetEventBinValues( std::vector<double> & vars )
{

   double logp = ( GetLeptonP() != 0 ? log10( GetLeptonP() ): 0 );

   vars.push_back( logp         );
   vars.push_back( GetCosineZ() );

}          

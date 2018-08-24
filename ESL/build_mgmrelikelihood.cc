#include <iostream>
#include <map>
#include <string>
#include <sstream>

#include "TFile.h"
#include "TTree.h"
#include "TChain.h"

#include "DataManager.h"
#include "CardReader.h"
#include "FileRecord.h"
#include "mgmreLikelihood.h"


int main( int argc, char * argv[] )
{

   if( argc < 5 ) 
   {
       std::cout << " -- Usage: " << std::endl;
       std::cout << " > " << argv[0] << " skx card cpu ncpus " << std::endl;
       std::cout << "    where: "  <<std::endl;
       std::cout << "           skx   = [ sk1, sk2, sk3, sk4 ] " << std::endl;
       std::cout << "           card  = Input card file " << std::endl;
       std::cout << "           cpu   = number of this cpu [0, ncpus -1] " << std::endl;
       std::cout << "           ncpus = total number of cpus used to process " << std::endl;
       exit( -1 ); 
   }
   else 
   {
      std::cout << " Command: " << std::endl;
      for( int i = 0 ; i < argc ; i ++ )
         std::cout << " " << argv[i] ; 

      std::cout << std::endl;
   }

   std::string mode  = "fcmc" ;
   std::string sk_era ( argv[1] );  // which detector geometry
   std::string card   ( argv[2] );
   int this_cpu = atoi( argv[3] );
   int tot_cpus = atoi( argv[4] );
   int seed;                        // where in the list of input files to start
   int blossom;                     // where to stop
   int skgen = -1;                  // sk generation {1,...4,...woa}


   if( sk_era.find("1") != std::string::npos ) skgen = 1; 
   if( sk_era.find("2") != std::string::npos ) skgen = 2; 
   if( sk_era.find("3") != std::string::npos ) skgen = 3; 
   if( sk_era.find("4") != std::string::npos ) skgen = 4; 

   skgen -= 1;
   if( skgen < 0 ) 
   { 
      std::cout << "Unrecognized SK Generation: " << sk_era << std::endl;
      exit(-1);
   }
   else 
      std::cout << sk_era << " corresponds to skgen = " << skgen <<std::endl;

   // load in all the data from the input card
   CardReader * MasterCard = new CardReader( card.c_str() );

   std::stringstream ss;
   std::string output_file;
   ss.str(""); ss << "output_mgmrelikelihood" ;
   MasterCard->GetKey( ss.str().c_str() , output_file );

   int kUseFiTQun = -1;
   MasterCard->GetKey( "use_fitqun",  kUseFiTQun );


   // Chain and load the files
   TChain * lchain = FileRecord::CheckAndBuildChain( MasterCard, mode  );
   int nEntries = (int) lchain->GetEntries();
   
   // Establish how many entries in the tree should be used 
   // per instance of the program (useful for batch mode)
   seed = this_cpu * int( nEntries / tot_cpus );
   blossom = ( seed+ int( nEntries / tot_cpus ) >= nEntries ? 
                     nEntries : seed + int( nEntries / tot_cpus ) );
   if ( this_cpu == tot_cpus - 1 ) blossom = nEntries;


   // Use DataManager to automagically read 
   // in and load all of the Branch information 
   DataManager * dm = new DataManager( lchain );
   dm->SetAllBranches( );

   mgmreLikelihood * llmgmre = new mgmreLikelihood( dm , skgen );

   
   // automatically add suffix to  
   // specified output file string in case 
   // we are processed in batch mode
   ss.str(""); 
   ss << output_file.substr(0, output_file.length() - 5 );
   if( tot_cpus > 1 ) 
     ss << "." << this_cpu << ".root";
   else 
     ss << ".root" ;

   std::cout << " Output will be written to: " << ss.str() << std::endl;
   llmgmre->SetOutputFile( ss.str().c_str() ); 

   // now we can loop over the input tree (or chain)
   std::cout << " Will process entries [" << seed << " , " << blossom << ") "	<< std::endl;

   llmgmre->ConstructLikelihood();

   // finish and close
   llmgmre->WriteHistograms();

   delete lchain;

   return 0;
}



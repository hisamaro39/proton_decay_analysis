#include <iostream>
#include <map>
#include <string>
#include <sstream>
#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "TSystem.h"
#include "tools/DataManager.h"
#include "tools/CardReader.h"
#include "src/OscNtupleManager.h"
#include "tools/TypeWrapper.h"
#include "core/FileRecord.h"

TChain * CheckAndBuildChain(  CardReader *, std::string );

int main( int argc, char * argv[] )
{
  if( argc != 8)
    {
      std::cout << " -- Usage: " << std::endl;
      std::cout << " > ./build_osc_ntuple skx mode card cpu ncpus (separation_setfile funcread_io func_file) " << std::endl;
      std::cout << "    where: "  <<std::endl;
      std::cout << "           skx   = [ sk1, sk2, sk3, sk4 ] " << std::endl;
      std::cout << "           input  = [ fcdt, fcmc, p_epi, p_mupi, p_eee, p_mumumu, p_emumu, p_muee ] " << std::endl;
      std::cout << "           mode  = [ p_epi, p_mupi, p_eee, p_mumumu, p_emumu, p_muee ] " << std::endl;
      std::cout << "           card  = Input card file " << std::endl;
      std::cout << "           cpu   = number of this cpu [0, ncpus -1] " << std::endl;
      std::cout << "           ncpus = total number of cpus used to process " << std::endl;
      std::cout << "           batch = using batch or not " << std::endl;
      std::cout << std::endl;
      exit( -1 );
    }   
   else 
   {
      std::cout << " Command: " << std::endl;
      for( int i = 0 ; i < argc ; i ++ )
         std::cout << " " << argv[i] ; 

      std::cout << std::endl;
   }

   std::string sk_era ( argv[1] );  // which detector geometry
   std::string input   ( argv[2] );  // fcdt, pcdt, fcmc, pcmc, ummc, 
   std::string mode   ( argv[3] );  // fcdt, pcdt, fcmc, pcmc, ummc, 
   std::string card   ( argv[4] );
   std::string this_cpu_str   ( argv[5] );
   int this_cpu = atoi( argv[5] );
   int tot_cpus = atoi( argv[6] );
   int use_batch = atoi( argv[7] );
   int seed;                        // where in the list of input files to start
   int blossom;                     // where to stop
   int skgen = -1;                  // sk generation {1,...4,...woa}


   if( sk_era.find("1") != std::string::npos ) skgen = 1; 
   if( sk_era.find("2") != std::string::npos ) skgen = 2; 
   if( sk_era.find("3") != std::string::npos ) skgen = 3; 
   if( sk_era.find("4") != std::string::npos ) skgen = 4; 
   if( skgen < 1 ) 
   { 
      std::cout << "Unrecognized SK Generation: " << sk_era << std::endl;
      exit(-1);
   }
   else 
      std::cout << sk_era << " corresponds to skgen = " << skgen <<std::endl;

   // load in all the data from the input card
   CardReader * MasterCard = new CardReader( card.c_str() );

   int kUseFiTQun = -1;
   MasterCard->GetKey( "use_fitqun",  kUseFiTQun );
   int kUseTauNN = -1;
   MasterCard->GetKey( "use_taunn", kUseTauNN);
   int kDebugMode = -1;
   MasterCard->GetKey( "debug_mode", kDebugMode);
   int kCheckBkg = -1;
   MasterCard->GetKey( "check_bkg", kCheckBkg);
   int kAllHist = -1;
   MasterCard->GetKey( "all_hist", kAllHist);
   int kCorrelatedDecay = -1;
   MasterCard->GetKey( "correlated_decay", kCorrelatedDecay);
   int kFermiMotion = -1;
   MasterCard->GetKey( "fermi_motion", kFermiMotion);
   int kOutsideSR = -1;
   MasterCard->GetKey( "outside_sr", kOutsideSR);
   int kCR = -1;
   MasterCard->GetKey( "use_cr", kCR);
   int kTotalBox = -1;
   MasterCard->GetKey( "only_total_box", kTotalBox);
   int kLiveTime = -1;
   MasterCard->GetKey( "use_live_time", kLiveTime);
   int kAllWeight = -1;
   MasterCard->GetKey( "use_all_weight", kAllWeight);
   int kSystNtag = -1;
   MasterCard->GetKey( "syst_ntag", kSystNtag);
   int kEnergyScale = -1;
   MasterCard->GetKey( "energy_scale", kEnergyScale);
   int kNonUni = -1;
   MasterCard->GetKey( "non_uniformity", kNonUni);
   int kPID = -1;
   MasterCard->GetKey( "pid", kPID);
   int kMakeNtuple = -1;
   MasterCard->GetKey( "make_ntuple", kMakeNtuple);
   std::string output_file,output_ntuple;
   /*if(use_batch) {
     if(kCheckBkg) output_file = "output_batch/" + input + "_" + mode + "/" + sk_era + "/file/" + input + "." + sk_era + "." + "mode_" + mode + "_check_bkg." + this_cpu_str + ".root";
     else output_file = "output_batch/" + input + "_" + mode + "/" + sk_era + "/file/" + input + "." + sk_era + "." + "mode_" + mode + "." + this_cpu_str + ".root";
   }*/
   //else {
   if(use_batch) output_file = "output_batch/" + input + "_" + mode + "/" + sk_era + "/file/" + input + "." + sk_era + "." + "mode_" + mode;
   else output_file = "output/" + input + "." + sk_era + "." + "mode_" + mode;
   if(use_batch) output_ntuple = "output_batch/" + input + "_" + mode + "/" + sk_era + "/file/" + input + "." + sk_era + "." + "mode_" + mode;
   else output_ntuple = "output/" + input + "." + sk_era + "." + "mode_" + mode;
   if(kCorrelatedDecay==0) output_file += "_cddown";
   if(kCorrelatedDecay==2) output_file += "_cdup";
   if(kFermiMotion) output_file += "_fermigas";
   if(kOutsideSR) output_file += "_outsideSR";
   if(kCR) output_file += "_CR";
   if(kCR) output_ntuple += "_CR";
   if(kTotalBox) output_ntuple += "_total_box";
   stringstream iteration;
   iteration << kSystNtag;
   if(kSystNtag) output_file += "_syst_ntag_it" + iteration.str();
   if(kEnergyScale==1) output_file += "_energydown_sk1";
   if(kEnergyScale==2) output_file += "_energyup_sk1";
   if(kEnergyScale==3) output_file += "_energydown_sk2";
   if(kEnergyScale==4) output_file += "_energyup_sk2";
   if(kEnergyScale==5) output_file += "_energydown_sk3";
   if(kEnergyScale==6) output_file += "_energyup_sk3";
   if(kEnergyScale==7) output_file += "_energydown_sk4";
   if(kEnergyScale==8) output_file += "_energyup_sk4";
   if(kNonUni==1) output_file += "_nonunidown_sk1";
   if(kNonUni==2) output_file += "_nonuniup_sk1";
   if(kNonUni==3) output_file += "_nonunidown_sk2";
   if(kNonUni==4) output_file += "_nonuniup_sk2";
   if(kNonUni==5) output_file += "_nonunidown_sk3";
   if(kNonUni==6) output_file += "_nonuniup_sk3";
   if(kNonUni==7) output_file += "_nonunidown_sk4";
   if(kNonUni==8) output_file += "_nonuniup_sk4";
   if(kPID==1) output_file += "_pid_sk1";
   if(kPID==2) output_file += "_pid_sk2";
   if(kPID==3) output_file += "_pid_sk3";
   if(kPID==4) output_file += "_pid_sk4";
   if(!kAllHist) output_file += "_validation";
   if(kCheckBkg) output_file += "_check_bkg";
   if(!kLiveTime) output_file += "_wo_livetime";
   if(!kAllWeight) output_file += "_wo_weight";
   if(use_batch) output_file += "." + this_cpu_str + ".root";
   else output_file += ".root";
   if(use_batch) output_ntuple += "_tree." + this_cpu_str + ".root";
   else output_ntuple += "_tree.root";
   if(kDebugMode || kMakeNtuple) output_file = "output/temp.root";
   //if(use_batch) output_ntuple = "output_batch/" + input + "_" + mode + "/" + sk_era + "/file/" + input + "." + sk_era + "." + "mode_" + mode + "_tree." + this_cpu_str + ".root";
   //else output_ntuple = "output/" + input + "." + sk_era + "." + "mode_" + mode + "_tree.root";


   std::cout << "output file is " << output_file << std::endl;
   std::cout << "output ntuple is " << output_ntuple << std::endl;

   float live_time_fcmc=-1, live_time_fcdt=-1;
   MasterCard->GetKey( "live_time_fcmc", live_time_fcmc);
   MasterCard->GetKey( "live_time_fcdt", live_time_fcdt);
   float weight_live_time = (input.find("fcmc")!=std::string::npos)? live_time_fcdt / live_time_fcmc : 1.;
   std::cout << "Live Time fcdt/fcmc/weight=" 
     << live_time_fcdt << "/" << live_time_fcmc << "/" << weight_live_time<< std::endl;

   // Chain and load the files
   std::cout << "Chain and load the files" << std::endl;
   TChain * lchain = FileRecord::CheckAndBuildChain( MasterCard, input  );
   int nEntries = (int) lchain->GetEntries();
   std::cout << "nEntries=" << nEntries << std::endl;

   //SKEventParser *Parser  =new SKEventParser();
   
   // Establish how many entries in the tree should be used 
   // per instance of the program (useful for batch mode)
   seed = this_cpu * int( nEntries / tot_cpus );
   blossom = ( seed+ int( nEntries / tot_cpus ) >= nEntries ? 
                     nEntries : seed + int( nEntries / tot_cpus ) );
   if ( this_cpu == tot_cpus - 1 ) blossom = nEntries;

   
   // Use DataManager to automagically read 
   // in and load all of the Branch information 
   std::cout << "Use DataManager" << std::endl;
   DataManager * dm = new DataManager( lchain );
   dm->SetAllBranches( );
   //Parser->LoadDataManager(dm);

   // Will perform the event classification
   //This uses the default card
   // OscNtupleManager * om  = new OscNtupleManager( dm , skgen ); 
   //This uses the specified card
   std::cout << "Define OscNtupleManager om" << std::endl;
   OscNtupleManager * om  = new OscNtupleManager( dm , MasterCard, skgen , output_file.c_str()); 
   
   std::cout << " mode       : " << mode << std::endl;
   std::cout << " kUseFiTQun : " << kUseFiTQun << std::endl;
   std::cout << " kUseTauNN : " << kUseTauNN << std::endl;
   std::cout << " kAllHist : " << kAllHist << std::endl;

   om->SetMode( mode  );
   om->SetInput( input  );
   if(kMakeNtuple) om->SetOutputNtuple( output_ntuple  );
   else om->SetOutputHist( output_file  );
   om->UseFiTQun(  (kUseFiTQun == 1 ? true : false ) );
   om->UseTauNN( (kUseTauNN == 1 ? true : false ) );
   om->DebugMode( (kDebugMode == 1 ? true : false ) );
   om->CheckBkg( (kCheckBkg == 1 ? true : false ) );
   om->AllHist( (kAllHist == 1 ? true : false ) );
   om->MakeNtuple( (kMakeNtuple == 1 ? true : false ) );
   om->SetInputTree( lchain );
   om->SetLiveTimeWeight(weight_live_time);
   om->CorrelatedDecay(kCorrelatedDecay );
   om->FermiMotion(kFermiMotion );
   om->OutsideSR(kOutsideSR );
   om->UseCR(kCR );
   om->TotalBox(kTotalBox );
   om->UseLiveTime(kLiveTime );
   om->UseAllWeight(kAllWeight );
   om->SystNtag(kSystNtag );
   om->EnergyScale(kEnergyScale );
   om->NonUni(kNonUni );
   om->PID(kPID );

   // automatically add suffix to  
   // specified output file string in case 
   // we are processed in batch mode
   std::stringstream ss;
   ss.str(""); 
   ss << output_file.substr(0, output_file.length() - 5 );
   if( tot_cpus > 1 ) 
     ss << "." << this_cpu << ".root";
   else 
     ss << ".root" ;

   std::cout << " Output will be written to: " << ss.str() << std::endl;

   // now we can loop over the input tree (or chain)
   std::cout << " Will process entries [" << seed << " , " << blossom << ") "	<< std::endl;
   om->Process( seed, blossom); 
   
   delete lchain;

   return 0;
}

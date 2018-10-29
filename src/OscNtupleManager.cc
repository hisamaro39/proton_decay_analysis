#include "src/OscNtupleManager.h"
#include <fstream>

#include "TLorentzVector.h"

#include "skcore/global.h"


using namespace std;


const std::string proton_type[2] = {"bound_proton","free_proton"};
const int range_match_gamma_mom[] = {0,200,300,400,500};
const int n_range_match_gamma_mom = sizeof(range_match_gamma_mom)/sizeof(int);
const int range_match_e_mom[] = {0,200,300,400,500};
const int n_range_match_e_mom = sizeof(range_match_e_mom)/sizeof(int);
const int range_match_mu_mom[] = {0,200,300,400,500};
const int n_range_match_mu_mom = sizeof(range_match_mu_mom)/sizeof(int);
const int range_momentum[] = {100,250,400,630,1000,2500,5000,10000,100000};
const int n_range_momentum = sizeof(range_momentum)/sizeof(int);
const float eff_e_ring[] = {0,0.06,0.15,0.34,0.65,0.75,0.86,0.92,0.95,0.97};//20GeV 


OscNtupleManager::OscNtupleManager( DataManager * dman , int skx , const char *filename) : m_hSvc(filename)

{
  dm     = dman;
  skgen  = skx-1;

  llm = new LikelihoodManager( skx );
  m_hSvc.addTrigger("");
  m_hSvc.addSystematics("");
  Initialize();
}


OscNtupleManager::OscNtupleManager( DataManager * dman , CardReader * cr , int skx , const char *filename) : m_hSvc(filename)

{
  dm     = dman;
  skgen  = skx-1;

  llm = new LikelihoodManager( cr, skx );

  m_hSvc.addTrigger("");
  m_hSvc.addSystematics("");

  Initialize();
}


void OscNtupleManager::Initialize() 
{

  os = new outputOscStructure;
  llm->SetDataManager( dm );

  nhitac_cut[0] =  10; 
  nhitac_cut[1] =  16; 
  nhitac_cut[2] =  16; 
  nhitac_cut[3] =  16; 

  potot_cut[0] =  3.0e3; 
  potot_cut[1] =  1.5e3; 
  potot_cut[2] =  3.0e3; 
  potot_cut[3] =  3.0e3; 

  llMGMRE       = -1e7 ;
  llMGMRENUE    = -1e7 ;
  llMGMRENUEBar = -1e7 ;

  kUseFiTQun    = false ;
  kUseTauNN     = false ; 
  kDebugMode     = false ; 
  kAllHist     = false ; 
  kMakeNtuple     = false ; 
  fsiweight     = 0 ;

  LoadWrappers();

  //m_fort->Initialize();

  m_fort = new fort_func();
  bNu = new BargerPropagator();

  //count event type 
  n_free_proton=0;
  //p -> e+ e+ e-
  n_eee=0;n_eeeg=0;n_eeep=0;n_eeen=0;n_eeegp=0;n_eeenp=0;n_eeegn=0;
  //p -> e+ pi0
  n_noint=0;n_abs=0;n_scat=0;n_charge=0;n_prod=0;n_match_e=0;n_match_0gamma=0;n_match_1gamma=0;n_match_2gamma=0;
  //p -> mu+ pi0
  n_match_mu=0;

  expected_3ring_events_electron=0.;
  expected_3ring_events_muon=0.;

  graph_point=0;graph_point_2=0;


}


OscNtupleManager::~OscNtupleManager( )
{
  delete os;
  delete llm;


  // delete llmgmre;
  // delete llpi0;
}

void OscNtupleManager::CreateBranch(){

  otree = new TTree("osc_tuple", "tree build for SK osc analyses"); 
  otree->Branch("ipnu"   , &o_ipnu    , "ipnu/I"    );
  otree->Branch("itype"   , &o_itype    , "itype/I"    );
  otree->Branch("pnu"   , &o_pnu    , "pnu/F"    );
  otree->Branch("dir"    , &o_dir     , "dir[3]/F"  );
  otree->Branch("amom"   , &o_amom    , "amom/F"    );
  otree->Branch("mode"   , &o_mode    , "mode/I"    );

}

void OscNtupleManager::CreateHist()
{

  //oscillation plot
  if(process_mode=="fcmc") {
    for(int i=1;i<15;i++){//interacion type
      m_hSvc.create1D(Form("nRing_type%d",i),"",10,0,10);
      m_hSvc.create1D(Form("nRing_type%d_osc",i),"",10,0,10);
      m_hSvc.create1D(Form("zenith_angle_type%d_osc",i),"",10,-1,1);
      m_hSvc.create1D(Form("momentum_log10_type%d",i),"",40,1,5);
      m_hSvc.create1D(Form("momentum_log10_type%d_osc",i),"",40,1,5);
      for(int m=0;m<n_range_momentum-1;m++){
        m_hSvc.create1D(Form("zenith_angle_type%d_mom%d_%d",i,range_momentum[m],range_momentum[m+1]),"",10,-1,1);
        m_hSvc.create1D(Form("zenith_angle_type%d_mom%d_%d_osc",i,range_momentum[m],range_momentum[m+1]),"",10,-1,1);
      }
    }
    return;
  }

  //validation plot
  if(!kAllHist){
    m_hSvc.createGraph("tgraph_test","");
    m_hSvc.create1D("true_mom_lepton","",40,0,800);
    m_hSvc.create1D("true_mom_lepton_match_ring","",40,0,800);
    m_hSvc.create1D("true_mom_lepton_match_ring_angle_elike","",40,0,800);
    m_hSvc.create1D("true_mom_lepton_match_ring_angle_mulike","",40,0,800);
    m_hSvc.create1D("true_mom_lepton_match_ring_charge_elike","",40,0,800);
    m_hSvc.create1D("true_mom_lepton_match_ring_charge_mulike","",40,0,800);
    m_hSvc.create1D("true_angle_lepton_and_lepton","",36,0,180);
    m_hSvc.create1D("true_angle_lepton_and_lepton_match_ring","",36,0,180);
    m_hSvc.create1D("residual_emom","",80,-0.4,0.4);
    m_hSvc.create1D("residual_mmom","",80,-0.4,0.4);
    m_hSvc.create1D("residual_total_mass","",80,-0.4,0.4);
    m_hSvc.create1D("residual_total_mom","",80,-0.8,0.8);
    m_hSvc.create1D("residual_total_gamma_mass","",80,-0.4,0.4);
    m_hSvc.create1D("residual_total_gamma_mom","",80,-0.8,0.8);
    m_hSvc.create1D("diff_opening_angle","",80,-10,10);
    m_hSvc.create1D("diff_opening_angle_angle_elike","",80,-10,10);
    m_hSvc.create1D("diff_opening_angle_charge_elike","",80,-10,10);
    m_hSvc.create1D("diff_opening_angle_angle_mulike","",80,-10,10);
    m_hSvc.create1D("diff_opening_angle_charge_mulike","",80,-10,10);
    m_hSvc.create1D("diff_opening_angle_ip2","",80,-10,10);
    m_hSvc.create1D("diff_opening_angle_ip3","",80,-10,10);
    m_hSvc.create1D("diff_vertex_r","",100,0,100);
    m_hSvc.create1D("opening_angle","",100,0,50);
    m_hSvc.create1D("expected_opening_angle","",100,0,50);
    m_hSvc.create1D("n_true_decayE","",5,0,5);
    m_hSvc.create1D("nDecayE_true_decayE_1","",5,0,5);
    m_hSvc.create1D("nDecayE_true_decayE_2","",5,0,5);
    m_hSvc.create1D("nDecayE_true_decayE_3","",5,0,5);
    m_hSvc.create2D("true_mom_expected_opening_angle","",100,0,800,100,0,50);
    m_hSvc.create2D("true_mom_opening_angle","",100,0,800,100,0,50);
    m_hSvc.create1D("true_min_mom_lepton","",40,0,800);
    m_hSvc.create1D("true_mid_mom_lepton","",40,0,800);
    m_hSvc.create1D("true_max_mom_lepton","",40,0,800);
    m_hSvc.create1D("true_mom_muon","",40,0,800);
    m_hSvc.create1D("true_mom_electron","",40,0,800);
    m_hSvc.create1D("true_mom_gamma","",40,0,800);
    for(int p=0;p<6;p++){
      m_hSvc.create1D(Form("diff_opening_angle_mom%d_%d",100*p,100+100*p),"",80,-10,10);
      m_hSvc.create1D(Form("diff_opening_angle_muon_mom%d_%d",100*p,100+100*p),"",80,-10,10);
      m_hSvc.create1D(Form("diff_opening_angle_electron_mom%d_%d",100*p,100+100*p),"",80,-10,10);
      m_hSvc.create1D(Form("diff_opening_angle_diff_vertex_r%d_%d",5*p,5+5*p),"",80,-10,10);
      m_hSvc.create1D(Form("residual_emom_mom%d_%d",100*p,100+100*p),"",80,-0.4,0.4);
      m_hSvc.create1D(Form("residual_mmom_mom%d_%d",100*p,100+100*p),"",80,-0.4,0.4);
      m_hSvc.create1D(Form("prob_angle_electron_mom%d_%d",100*p,100+100*p),"",100,-50,50);
      m_hSvc.create1D(Form("prob_angle_muon_mom%d_%d",100*p,100+100*p),"",100,-5,5);
    }

    for(int r=2;r<4;r++){
      m_hSvc.create1D(Form("true_max_mom_lepton_nring%d",r),"",40,0,800);
      m_hSvc.create1D(Form("true_mid_mom_lepton_nring%d",r),"",40,0,800);
      m_hSvc.create1D(Form("true_min_mom_lepton_nring%d",r),"",40,0,800);
      m_hSvc.create1D(Form("true_min_angle_lepton_lepton_nring%d",r),"",36,0,180);
      m_hSvc.create1D(Form("true_angle_min_mid_lepton_nring%d",r),"",36,0,180);
      m_hSvc.create1D(Form("true_angle_min_max_lepton_nring%d",r),"",36,0,180);
      m_hSvc.create1D(Form("true_angle_mid_max_lepton_nring%d",r),"",36,0,180);
      m_hSvc.create1D(Form("true_mom_lepton_nring%d",r),"",40,0,800);
      m_hSvc.create1D(Form("true_mom_muon_nring%d",r),"",40,0,800);
      m_hSvc.create1D(Form("true_mom_lepton_match_ring_nring%d",r),"",40,0,800);
      m_hSvc.create1D(Form("true_mom_muon_match_ring_nring%d",r),"",40,0,800);
      m_hSvc.create1D(Form("true_angle_lepton_and_ring_nring%d",r),"",36,0,180);
      m_hSvc.create1D(Form("true_mom_lepton_match_ring_angle_elike_nring%d",r),"",40,0,800);
      m_hSvc.create1D(Form("true_mom_lepton_match_ring_angle_mulike_nring%d",r),"",40,0,800);
      m_hSvc.create1D(Form("true_mom_muon_match_ring_angle_mulike_nring%d",r),"",40,0,800);
      for(int m1=0;m1<2;m1++){
        m_hSvc.create1D(Form("residual_total_mass_nring%d_mulike%d",r,m1),"",80,-0.4,0.4);
        m_hSvc.create1D(Form("residual_total_mom_nring%d_mulike%d",r,m1),"",80,-0.4,0.4);
        m_hSvc.create1D(Form("residual_total_gamma_mass_nring%d_mulike%d",r,m1),"",80,-0.4,0.4);
        m_hSvc.create1D(Form("residual_total_gamma_mom_nring%d_mulike%d",r,m1),"",80,-0.4,0.4);
        m_hSvc.create1D(Form("true_mom_muon_nring%d_mulike%d",r,m1),"",40,0,800);
        m_hSvc.create1D(Form("diff_vertex_r_nring%d_mulike%d",r,m1),"",100,0,100);
        m_hSvc.create1D(Form("nMulikeRing_angle_nring%d_trueDecayE%d",r,m1),"",6,0,6);
      }
    }
    for(int f=0;f<2;f++){//free proton or not
      m_hSvc.create1D(Form("total_true_mass_fp%d",f),"",125,0,1250);
      m_hSvc.create1D(Form("total_true_mom_fp%d",f),"",100,0,1000);
      m_hSvc.create2D(Form("total_true_mass_true_mom_fp%d",f),"",125,0,1250,100,0,1000);
      m_hSvc.create1D(Form("total_gamma_true_mass_fp%d",f),"",125,0,1250);
      m_hSvc.create1D(Form("total_gamma_true_mom_fp%d",f),"",100,0,1000);
      m_hSvc.create2D(Form("total_true_mass_1st_lepton_mom_fp%d",f),"",125,0,1250,100,0,800);
      m_hSvc.create2D(Form("total_true_mom_1st_lepton_mom_fp%d",f),"",100,0,1000,100,0,800);
      m_hSvc.create1D(Form("true_mom_1st_lepton_fp%d",f),"",40,0,800);
      m_hSvc.create1D(Form("true_mom_2nd_lepton_fp%d",f),"",40,0,800);
      m_hSvc.create1D(Form("true_mom_3rd_lepton_fp%d",f),"",40,0,800);
      m_hSvc.create1D(Form("true_energy_1st_lepton_fp%d",f),"",40,0,800);
      m_hSvc.create1D(Form("true_energy_2nd_lepton_fp%d",f),"",40,0,800);
      m_hSvc.create1D(Form("true_energy_3rd_lepton_fp%d",f),"",40,0,800);
      m_hSvc.create1D(Form("true_mass_1st_2nd_fp%d",f),"",50,0,1250);
      m_hSvc.create1D(Form("true_mass_2nd_3rd_fp%d",f),"",50,0,1250);
      m_hSvc.create1D(Form("true_mass_3rd_1st_fp%d",f),"",50,0,1250);
      m_hSvc.createGraph(Form("true_mass_1st_2nd_vs_2nd_3rd_fp%d",f),"");

    }
  }

  if(kAllHist){

    int r_max=3,mu_max=3,p_max=3;
    for(int c=0;c<10;c++){//cut type
      for(int r=0;r<r_max;r++){//cut pattern of nring
        for(int mu=0;mu<mu_max;mu++){//cut pattern of mu-like ring
          for(int p=0;p<p_max;p++){//cut pattern of michel electron
            if(kCheckBkg){//default is fcmc
              for(int im=1;im<93;im++){
                m_hSvc.create1D(Form("distance_to_wall_cut%d_nring%d_mulike%d_michel%d_mode%d",c,r,mu,p,im),"",100,0,1000);
                m_hSvc.create1D(Form("visible_energy_cut%d_nring%d_mulike%d_michel%d_mode%d",c,r,mu,p,im),"",100,0,500);
                m_hSvc.create1D(Form("nhit_OD_cut%d_nring%d_mulike%d_michel%d_mode%d",c,r,mu,p,im),"",100,0,100);
                m_hSvc.create1D(Form("nRing_cut%d_nring%d_mulike%d_michel%d_mode%d",c,r,mu,p,im),"",10,0,10);
                m_hSvc.create1D(Form("n_michel_electron_cut%d_nring%d_mulike%d_michel%d_mode%d",c,r,mu,p,im),"",5,0,5);
                for(int e=1;e<4;e++) m_hSvc.create1D(Form("mass_pi0_reco_elike%d_cut%d_nring%d_mulike%d_michel%d_mode%d",e,c,r,mu,p,im),"",50,0,500);
                m_hSvc.create1D(Form("ntag_multiplicity_cut%d_nring%d_mulike%d_michel%d_mode%d",c,r,mu,p,im),"",10,0,10);
                m_hSvc.create1D(Form("mass_proton_reco_cut%d_nring%d_mulike%d_michel%d_mode%d",c,r,mu,p,im),"",125,0,1250);
                m_hSvc.create1D(Form("mom_proton_reco_cut%d_nring%d_mulike%d_michel%d_mode%d",c,r,mu,p,im),"",100,0,1000);
              }
            }
            m_hSvc.create1D(Form("interaction_mode_cut%d_nring%d_mulike%d_michel%d",c,r,mu,p),"",92,1,93);
            m_hSvc.create1D(Form("distance_to_wall_cut%d_nring%d_mulike%d_michel%d",c,r,mu,p),"",100,0,1000);
            m_hSvc.create1D(Form("visible_energy_cut%d_nring%d_mulike%d_michel%d",c,r,mu,p),"",100,0,500);
            m_hSvc.create1D(Form("nhit_OD_cut%d_nring%d_mulike%d_michel%d",c,r,mu,p),"",100,0,100);
            m_hSvc.create1D(Form("nRing_cut%d_nring%d_mulike%d_michel%d",c,r,mu,p),"",10,0,10);
            m_hSvc.create1D(Form("n_michel_electron_cut%d_nring%d_mulike%d_michel%d",c,r,mu,p),"",5,0,5);
            m_hSvc.create1D(Form("ntag_multiplicity_cut%d_nring%d_mulike%d_michel%d",c,r,mu,p),"",10,0,10);
            m_hSvc.create1D(Form("mass_proton_reco_cut%d_nring%d_mulike%d_michel%d",c,r,mu,p),"",125,0,1250);
            m_hSvc.create1D(Form("mom_proton_reco_cut%d_nring%d_mulike%d_michel%d",c,r,mu,p),"",100,0,1000);
            m_hSvc.createGraph(Form("mass_mom_proton_reco_cut%d_nring%d_mulike%d_michel%d",c,r,mu,p),"");
            m_hSvc.createGraph(Form("all_ring_mass_mom_proton_reco_cut%d_nring%d_mulike%d_michel%d",c,r,mu,p),"");
            m_hSvc.createGraph(Form("all_mulike_mass_mom_proton_reco_cut%d_nring%d_mulike%d_michel%d",c,r,mu,p),"");
            for(int n=2;n<4;n++){
              m_hSvc.create1D(Form("nElikeRing_angle_nring%d_cut%d_nring%d_mulike%d_michel%d",n,c,r,mu,p),"",6,0,6);
              m_hSvc.create1D(Form("nElikeRing_pattern_nring%d_cut%d_nring%d_mulike%d_michel%d",n,c,r,mu,p),"",6,0,6);
              m_hSvc.create1D(Form("nMulikeRing_angle_nring%d_cut%d_nring%d_mulike%d_michel%d",n,c,r,mu,p),"",6,0,6);
              m_hSvc.create1D(Form("nMulikeRing_pattern_nring%d_cut%d_nring%d_mulike%d_michel%d",n,c,r,mu,p),"",6,0,6);
            }

            if(process_input=="fcmc" || process_input=="fcdt") continue;

            for(int f=0;f<2;f++){
              m_hSvc.create1D(Form("interaction_mode_cut%d_nring%d_mulike%d_michel%d_fp%d",c,r,mu,p,f),"",92,1,93);
              m_hSvc.create1D(Form("distance_to_wall_cut%d_nring%d_mulike%d_michel%d_fp%d",c,r,mu,p,f),"",100,0,1000);
              m_hSvc.create1D(Form("visible_energy_cut%d_nring%d_mulike%d_michel%d_fp%d",c,r,mu,p,f),"",100,0,500);
              m_hSvc.create1D(Form("nhit_OD_cut%d_nring%d_mulike%d_michel%d_fp%d",c,r,mu,p,f),"",100,0,100);
              m_hSvc.create1D(Form("nRing_cut%d_nring%d_mulike%d_michel%d_fp%d",c,r,mu,p,f),"",10,0,10);
              m_hSvc.create1D(Form("n_michel_electron_cut%d_nring%d_mulike%d_michel%d_fp%d",c,r,mu,p,f),"",5,0,5);
              m_hSvc.create1D(Form("ntag_multiplicity_cut%d_nring%d_mulike%d_michel%d_fp%d",c,r,mu,p,f),"",10,0,10);
              m_hSvc.create1D(Form("mass_proton_reco_cut%d_nring%d_mulike%d_michel%d_fp%d",c,r,mu,p,f),"",125,0,1250);
              m_hSvc.create1D(Form("mom_proton_reco_cut%d_nring%d_mulike%d_michel%d_fp%d",c,r,mu,p,f),"",100,0,1000);
              m_hSvc.createGraph(Form("mass_mom_proton_reco_cut%d_nring%d_mulike%d_michel%d_fp%d",c,r,mu,p,f),"");
              m_hSvc.createGraph(Form("all_ring_mass_mom_proton_reco_cut%d_nring%d_mulike%d_michel%d_fp%d",c,r,mu,p,f),"");
              m_hSvc.createGraph(Form("all_mulike_mass_mom_proton_reco_cut%d_nring%d_mulike%d_michel%d_fp%d",c,r,mu,p,f),"");
              for(int n=2;n<4;n++){
                m_hSvc.create1D(Form("nElikeRing_angle_nring%d_cut%d_nring%d_mulike%d_michel%d_fp%d",n,c,r,mu,p,f),"",6,0,6);
                m_hSvc.create1D(Form("nElikeRing_pattern_nring%d_cut%d_nring%d_mulike%d_michel%d_fp%d",n,c,r,mu,p,f),"",6,0,6);
                m_hSvc.create1D(Form("nMulikeRing_angle_nring%d_cut%d_nring%d_mulike%d_michel%d_fp%d",n,c,r,mu,p,f),"",6,0,6);
                m_hSvc.create1D(Form("nMulikeRing_pattern_nring%d_cut%d_nring%d_mulike%d_michel%d_fp%d",n,c,r,mu,p,f),"",6,0,6);
                m_hSvc.create1D(Form("mass_all_ring_reco_nring%d_cut%d_nring%d_mulike%d_michel%d_fp%d",n,c,r,mu,p,f),"",125,0,1250);
                m_hSvc.create1D(Form("mom_all_ring_reco_nring%d_cut%d_nring%d_mulike%d_michel%d_fp%d",n,c,r,mu,p,f),"",125,0,1250);
                m_hSvc.create1D(Form("mass_all_mulike_reco_nring%d_cut%d_nring%d_mulike%d_michel%d_fp%d",n,c,r,mu,p,f),"",125,0,1250);
                m_hSvc.create1D(Form("mom_all_mulike_reco_nring%d_cut%d_nring%d_mulike%d_michel%d_fp%d",n,c,r,mu,p,f),"",125,0,1250);
              }
            }
          }
        }
      }
    }

    for(int r=0;r<r_max;r++){
      for(int mu=0;mu<mu_max;mu++){
        for(int m=0;m<p_max;m++){
          m_hSvc.create1D(Form("cut_flow_nring%d_mulike%d_michel%d",r,mu,m),"",11,0,11);
          for(int f=0;f<2;f++){
            m_hSvc.create1D(Form("cut_flow_nring%d_mulike%d_michel%d_fp%d",r,mu,m,f),"",11,0,11);
          }
        }
      }
    }

  }//all hist

}

void OscNtupleManager::LoadWrappers()
{

  dm->Get( "pnu"  , pnu    );
  dm->Get( "ndcy"  , ndcy    );
  dm->Get( "ngate"  , ngate    );
  dm->Get( "nbye"  , nbye    );
  dm->Get("evis"  , evis   );
  dm->Get("wall"  , wall   );
  dm->Get("nring" , nring  );
  dm->Get("numnu" , numnu  );
  dm->Get("npar" , npar  );
  dm->Get("npar2" , npar2  );
  dm->Get("nhitac", nhitac );
  dm->Get("potot" , potot  );
  dm->Get("ip"    , ip     );
  dm->Get("ipv"    , ipv     );
  dm->Get("ipv2"    , ipv2     );
  dm->Get("iorg"    , iorg     );
  dm->Get("pos" , pos  );
  dm->Get("posv" , posv  );
  dm->Get("ang" , ang  );
  dm->Get("ange" , ange  );
  dm->Get("angm" , angm  );
  dm->Get("amom" , amom  );
  dm->Get("amome" , amome  );
  dm->Get("amomm" , amomm  );
  dm->Get("prmslg", prmslg );
  dm->Get("probms", probms );
  dm->Get("dir"   , dir    );
  dm->Get("dirv"   , dirv    );
  dm->Get("dirv2"   , dirv2    );
  dm->Get("msdir" , msdir  );
  dm->Get("catpc" , catpc  );
  dm->Get("catpc_qrat_corrected", catpc_qrat_corrected);

  dm->Get("dirtotmue", dirtotmue  );
  dm->Get("etotmue"  , etotmue    );
  dm->Get("dirnu"    , dirnu      );
  dm->Get("ipnu"     , ipnu       );
  dm->Get("mode"     , mode       );

  dm->Get("flxb03"   , flxb03     );

  dm->Get("flxh06"   , flxh06     );
  dm->Get("flxh11"   , flxh11     );

  dm->Get("oscwgt"   , oscwgt     );

  // needed by UpMu
  dm->Get("wallv"    , wallv     );
  dm->Get("wallv2"    , wallv2     );
  dm->Get("pmomv"    , pmomv     );
  dm->Get("pmomv2"    , pmomv2     );
  dm->Get("Fit_dir"  , fit_dir   );
  dm->Get("Fit_pid"  , fit_pid   );
  dm->Get("Fit_len"  , fit_len   );
  dm->Get("Fit_mom"  , fit_mom   );
  dm->Get("Sh_id"    , sh_id     );
  dm->Get("Sh_delta" , sh_delta  );
  dm->Get("Sh_mean"  , sh_mean   );
  dm->Get("Sh_meanq" , sh_meanq  );
  dm->Get("Sh_chi1p" , sh_chi1p  );

  // I think um_* may not be filled 
  // in the ntuples for sk1,3 and 
  // only in sk2...?  -rvw
  dm->Get("Um_ehit8m", um_ehit8m );
  dm->Get("Um_ohit8m", um_ohit8m );


  /// fiTQun related
  dm->Get("fqpi0nll"  , fqpi0nll  );
  dm->Get("fq1rnll"   , fq1rnll   );
  dm->Get("fqpi0mass" , fqpi0mass );
  dm->Get("fq4rnll"   , fq4rnll   );

  //Tau_NN realted
  dm->Get("NN_selected", NN_selected);
  dm->Get("NN_output", NN_output);
  fsiweight = dm->GetTArrayF( "weights" );

  dm->Get( "ntag_nn"  , ntag_nn    );
  dm->Get( "ntag_mctruth_nn"  , ntag_mctruth_nn    );

  //variable from friend file
  dm->Get( "ErmsHax"  , ErmsHax    );
  dm->Get( "nEAveHax"  , nEAveHax    );


  /////////////////////////////////////
  dm->EnableUsedListOnly();

}

void OscNtupleManager::Process(int seed, int blossom )
{
  cout << "kAllHist=" << kAllHist << endl;

  if (kMakeNtuple) CreateBranch();
  else CreateHist();

  //define values in ntuple

  int entries = _ltree->GetEntries();
  int startEntry = seed;
  int endEntry = ( blossom < 0 ? entries : blossom ); 
  //int endEntry = 100;//temporal 

  bool fill;

  int Good = 0;
  int Bad = 0;


  for( int i = startEntry ; i < endEntry ; i ++ )
  {
    if(kDebugMode) std::cout << "Entry is " << i << std::endl;
    else if(i%10000==0) std::cout << "Entry is " << i << std::endl;

    _ltree->GetEntry(i);

    ZeroStructure();

    if(kMakeNtuple){
      FillNtuple();
      continue;
    }

    weight = 1.;
    mc_weight = 1.;
    osc_weight = 1.;
    true_mode=-1;
    if(process_input=="fcmc" || process_input=="fcdt"){
      bool sc1 = SetEventType(process_input,interaction_type);
      sc1 = SetTrueMode(true_mode); 
      if(process_input=="fcmc") {
        weight = live_time_weight;
        mc_weight = GetMCweight(); 
        //osc_weight = Oscillate();//3D calculation
        osc_weight = oscwgt(0);
        //weight *= weight*osc_weight;//osc_weight include mc_weight
        weight *= mc_weight;
      }
      if(process_mode=="fcmc"){
        MakeOscillationPlot();
        continue;
      }
    }

    if(kDebugMode) {
      std::cout << "live time weight is " << live_time_weight << std::endl;
      std::cout << "weight mc is " << mc_weight << std::endl;
      std::cout << "weight osc is " << osc_weight << std::endl;
      std::cout << "weight=" << weight << std::endl;
      std::cout << "mode=" << mode(0) << std::endl;
      std::cout << "npar/npar2=" << npar(0) << "/" << npar2(0) << std::endl;
      std::cout << "true_mode=" << true_mode << std::endl;
      for(int m=0;m<npar(0);m++){
        cout << "particle_" << m << " pid=" << ipv(m) << endl;
      }
      for(int t=0;t<npar2(0);t++){
        cout << "daughter particle from " << iorg(t) << " pid=" << ipv2(t) << endl;
      }
    }

    if(process_mode=="p_epi") Process_pepi();
    if(process_mode=="p_mupi") Process_pmupi();
    if(process_mode=="p_eee") Process_peee();
    if(process_mode=="p_mumumu") Process_pmumumu();
    if(process_mode=="p_emumu") Process_pemumu();
    if(process_mode=="p_muee") Process_pmuee();
    if(process_mode=="p_eemu") Process_peemu();
  }

  std::cout << " Classification: " << std::endl
    << "   Written: " << Good << std::endl
    << "   Skipped: " << Bad << std::endl
    << std::endl;

  cout << "expected 3ring events electron/muon=" 
    << expected_3ring_events_electron << "/" << expected_3ring_events_muon << endl;

  if(kMakeNtuple){
    TFile *file = new TFile(output_ntuple.c_str(),"recreate");
    otree->Write();
    file->Close();
  }

  else m_hSvc.WriteOutput();

}


void OscNtupleManager::Process_pepi(){

  //apply selection here

  if(!kDebugMode && wallv(0)<200) return;

  pass_cut[0][0]=true;
  if ( evis(0) > 30.0 && wall(0) > 200.0 && nhitac(0) <  nhitac_cut[skgen] ){//FC & FV
    pass_cut[1][0]=true;
  }

  //# of cherenkov ring
  if(nring(0)==2) pass_cut[2][0]=true;
  if(nring(0)==3) pass_cut[2][1]=true;
  //if(nring(0)==4) pass_cut[2][2]=true;
  if(nring(0)==2 || nring(0)==3) pass_cut[2][2]=true;

  //# of mu-like ring
  if(n_mulike_pattern==0) pass_cut[3][0]=true;
  if(n_mulike_pattern==1) pass_cut[3][1]=true;
  if(n_mulike_pattern==2) pass_cut[3][2]=true;
  if(n_mulike_pattern==3) pass_cut[3][3]=true;

  // michel (decay) electron cut
  if(nDecayE==0) pass_cut[4][0]=true;
  if(nDecayE==1) pass_cut[4][1]=true;
  if(nDecayE==2) pass_cut[4][2]=true;
  if(nDecayE==3) pass_cut[4][3]=true;

  vector<TLorentzVector> egamma_cand;
  egamma_cand.clear();
  for(int r=0;r<nRing;r++)
    if( ( prmslg(r,1) - prmslg(r,2) ) < 0 ) egamma_cand.push_back(GetTLorentzVectorRing(r,0));

  int gamma1_id=-1,gamma2_id=-1;
  float min_diff=99999999;
  if(egamma_cand.size()>1){
    for(unsigned int r1=0;r1<egamma_cand.size()-1;r1++){
      for(unsigned int r2=r1+1;r2<egamma_cand.size();r2++){
        float mass_pi0_reco = (egamma_cand[r1] + egamma_cand[r2]).M();
        float diff = fabs(mass_pi0_reco-135);
        if(diff < min_diff){
          min_diff = diff;
          closest_mass_pi0_reco = mass_pi0_reco;
          gamma1_id = r1;
          gamma2_id = r2;
        }
      }
    }
  }

  if(nRing==3 && closest_mass_pi0_reco>85 && closest_mass_pi0_reco<185) pass_cut[5][0]=true;
  if(nRing==2) pass_cut[5][0]=true;


  TLorentzVector total_vec;
  if(n_elike_pattern==3) total_vec = egamma_cand[0] + egamma_cand[1] + egamma_cand[2];
  if(n_elike_pattern==2) total_vec = egamma_cand[0] + egamma_cand[1];
  total_mass = total_vec.M();
  total_mom = total_vec.P();
  if(total_mass>800 && total_mass<1050 && total_mom<100){//total mass & low momentum
    pass_cut[6][0]=true;
  }
  if(total_mass>800 && total_mass<1050 && total_mom>100 && total_mom<250){//total mass & high momentum
    pass_cut[7][0]=true;
  }

  if(nNeutron==0) pass_cut[8][0]=true;

  TLorentzVector all_ring_vec;
  for(int r=0;r<nRing;r++){
    TLorentzVector temp_vec = 
    ( ( prmslg(r,1) - prmslg(r,2) ) < 0 ) ? GetTLorentzVectorRing(r,0) : GetTLorentzVectorRing(r,2);
    all_ring_vec += temp_vec;
  }
  all_ring_mass = all_ring_vec.M();
  all_ring_mom = all_ring_vec.P();

  //MakeCutFlowValidate();
  if(kAllHist) MakeCutFlow();
  else MakeValidationPlot();

}

void OscNtupleManager::Process_pmupi(){

  //apply selection here

  if(!kDebugMode && wallv(0)<200) return;

  pass_cut[0][0]=true;
  if ( evis(0) > 30.0 && wall(0) > 200.0 && nhitac(0) <  nhitac_cut[skgen] ){//FC & FV
    pass_cut[1][0]=true;
  }

  //# of cherenkov ring
  if(nring(0)==2) pass_cut[2][0]=true;
  if(nring(0)==3) pass_cut[2][1]=true;
  //if(nring(0)==4) pass_cut[2][2]=true;
  if(nring(0)==2 || nring(0)==3) pass_cut[2][2]=true;

  //# of mu-like ring
  if(n_mulike_pattern==0) pass_cut[3][0]=true;
  if(n_mulike_pattern==1) pass_cut[3][1]=true;
  if(n_mulike_pattern==2) pass_cut[3][2]=true;
  if(n_mulike_pattern==3) pass_cut[3][3]=true;

  // michel (decay) electron cut
  if(nDecayE==0) pass_cut[4][0]=true;
  if(nDecayE==1) pass_cut[4][1]=true;
  if(nDecayE==2) pass_cut[4][2]=true;
  if(nDecayE==3) pass_cut[4][3]=true;

  vector<TLorentzVector> gamma_cand;
  TLorentzVector mu_cand;
  gamma_cand.clear();
  for(int r=0;r<nRing;r++){
    if ( ( prmslg(r,1) - prmslg(r,2) ) < 0 ) gamma_cand.push_back(GetTLorentzVectorRing(r,0));
    else mu_cand = GetTLorentzVectorRing(r,2);
  }

  int gamma1_id=-1,gamma2_id=-1;
  float min_diff=99999999;
  if(gamma_cand.size()>1){
    for(unsigned int r1=0;r1<gamma_cand.size()-1;r1++){
      for(unsigned int r2=r1+1;r2<gamma_cand.size();r2++){
        float mass_pi0_reco = (gamma_cand[r1] + gamma_cand[r2]).M();
        float diff = fabs(mass_pi0_reco-135);
        if(diff < min_diff){
          min_diff = diff;
          closest_mass_pi0_reco = mass_pi0_reco;
          gamma1_id = r1;
          gamma2_id = r2;
        }
      }
    }
  }

  if(nRing==3 && closest_mass_pi0_reco>85 && closest_mass_pi0_reco<185) pass_cut[5][0]=true;
  if(nRing==2) pass_cut[5][0]=true;


  TLorentzVector total_vec;
  if(n_elike_pattern==2) total_vec = gamma_cand[0] + gamma_cand[1] + mu_cand;
  if(n_elike_pattern==1) total_vec = gamma_cand[0] + mu_cand;
  total_mass = total_vec.M();
  total_mom = total_vec.P();
  if(total_mass>800 && total_mass<1050 && total_mom<100){//total mass & low momentum
    pass_cut[6][0]=true;
  }
  if(total_mass>800 && total_mass<1050 && total_mom>100 && total_mom<250){//total mass & high momentum
    pass_cut[7][0]=true;
  }

  if(nNeutron==0) pass_cut[8][0]=true;

  TLorentzVector all_ring_vec;
  for(int r=0;r<nRing;r++){
    TLorentzVector temp_vec = 
    ( ( prmslg(r,1) - prmslg(r,2) ) < 0 ) ? GetTLorentzVectorRing(r,0) : GetTLorentzVectorRing(r,2);
    all_ring_vec += temp_vec;
  }
  all_ring_mass = all_ring_vec.M();
  all_ring_mom = all_ring_vec.P();

  //MakeCutFlowValidate();
  if(kAllHist) MakeCutFlow();
  else MakeValidationPlot();

}

void OscNtupleManager::Process_peee(){

  //apply selection here

  if(!kDebugMode && wallv(0)<200) return;

  pass_cut[0][0]=true;
  if ( evis(0) > 30.0 && wall(0) > 200.0 && nhitac(0) <  nhitac_cut[skgen] ){//FC & FV
    pass_cut[1][0]=true;
  }

  //# of cherenkov ring
  if(nring(0)==2) pass_cut[2][0]=true;
  if(nring(0)==3) pass_cut[2][1]=true;
  //if(nring(0)==4) pass_cut[2][2]=true;
  if(nring(0)==2 || nring(0)==3) pass_cut[2][2]=true;

  //# of mu-like ring
  if(n_mulike_angle==0) pass_cut[3][0]=true;
  if(n_mulike_angle==1) pass_cut[3][1]=true;
  if(n_mulike_angle==2) pass_cut[3][2]=true;
  if(n_mulike_angle==3) pass_cut[3][3]=true;

  // michel (decay) electron cut
  if(nDecayE==0) pass_cut[4][0]=true;
  if(nDecayE==1) pass_cut[4][1]=true;
  if(nDecayE==2) pass_cut[4][2]=true;
  if(nDecayE==3) pass_cut[4][3]=true;

  vector<TLorentzVector> e_cand;
  e_cand.clear();
  for(int r=0;r<nRing;r++)
    if( ( sqrt(fabs(probms(r,1)))-sqrt(fabs(probms(r,2))) ) < 0) e_cand.push_back(GetTLorentzVectorRing(r,1));

  pass_cut[5][0] = true;//no pi0 selection


  TLorentzVector total_vec;
  if(n_elike_angle==3) total_vec = e_cand[0] + e_cand[1] + e_cand[2];
  if(n_elike_angle==2) total_vec = e_cand[0] + e_cand[1];
  total_mass = total_vec.M();
  total_mom = total_vec.P();
  if(total_mass>800 && total_mass<1050 && total_mom<100){//total mass & low momentum
    pass_cut[6][0]=true;
  }
  if(total_mass>800 && total_mass<1050 && total_mom>100 && total_mom<250){//total mass & high momentum
    pass_cut[7][0]=true;
  }

  if(nNeutron==0) pass_cut[8][0]=true;

  TLorentzVector all_ring_vec;
  for(int r=0;r<nRing;r++){
    TLorentzVector temp_vec = 
    ( ( sqrt(fabs(probms(r,1)))-sqrt(fabs(probms(r,2))) ) < 0) ? GetTLorentzVectorRing(r,0) : GetTLorentzVectorRing(r,2);
    all_ring_vec += temp_vec;
  }
  all_ring_mass = all_ring_vec.M();
  all_ring_mom = all_ring_vec.P();

  if(kAllHist) MakeCutFlow();
  else MakeValidationPlot();

}

void OscNtupleManager::Process_pmumumu(){

  //apply selection here

  if(!kDebugMode && wallv(0)<200) return;

  pass_cut[0][0]=true;
  if ( evis(0) > 30.0 && wall(0) > 200.0 && nhitac(0) <  nhitac_cut[skgen] ){//FC & FV
    pass_cut[1][0]=true;
  }

  if(nring(0)==2) pass_cut[2][0]=true;
  if(nring(0)==3) pass_cut[2][1]=true;
  if(nring(0)==4) pass_cut[2][2]=true;

  if(n_mulike_angle==0) pass_cut[3][0]=true;
  if(n_mulike_angle==1) pass_cut[3][1]=true;
  if(n_mulike_angle==2) pass_cut[3][2]=true;
  if(n_mulike_angle==3) pass_cut[3][3]=true;

  // michel (decay) electron cut
  if(nDecayE==0) pass_cut[4][0]=true;
  if(nDecayE==1) pass_cut[4][1]=true;
  if(nDecayE==2) pass_cut[4][2]=true;
  if(nDecayE==3) pass_cut[4][3]=true;

  vector<TLorentzVector> mu_cand;
  mu_cand.clear();
  for(int r=0;r<nRing;r++){
    if( ( sqrt(fabs(probms(r,1)))-sqrt(fabs(probms(r,2))) ) > 0) mu_cand.push_back(GetTLorentzVectorRing(r,2));
  }

  pass_cut[5][0]=true;//no pi0 mass cut


  TLorentzVector total_vec;
  if(n_mulike_angle==3) total_vec = mu_cand[0] + mu_cand[1] + mu_cand[2];
  if(n_mulike_angle==2) total_vec = mu_cand[0] + mu_cand[1];
  total_mass = total_vec.M();
  total_mom = total_vec.P();
  if(total_mass>800 && total_mass<1050 && total_mom<100){//total mass & low momentum
    pass_cut[6][0]=true;
  }
  if(total_mass>800 && total_mass<1050 && total_mom>100 && total_mom<250){//total mass & high momentum
    pass_cut[7][0]=true;
  }

  if(nNeutron==0) pass_cut[8][0]=true;

  TLorentzVector all_ring_vec,all_mulike_vec;
  for(int r=0;r<nRing;r++){
    TLorentzVector temp_vec = 
      ( ( sqrt(fabs(probms(r,1)))-sqrt(fabs(probms(r,2))) ) < 0) ? GetTLorentzVectorRing(r,0) : GetTLorentzVectorRing(r,2);
    all_ring_vec += temp_vec;
    all_mulike_vec += GetTLorentzVectorRing(r,2);
  }
  all_ring_mass = all_ring_vec.M();
  all_ring_mom = all_ring_vec.P();
  all_mulike_mass = all_mulike_vec.M();
  all_mulike_mom = all_mulike_vec.P();

  if(all_mulike_mass>800 && all_mulike_mass<1050 && all_mulike_mom<100){//all_mulike mass & low momentum
    pass_cut[9][0]=true;
  }
  if(all_mulike_mass>800 && all_mulike_mass<1050 && all_mulike_mom>100 && all_mulike_mom<250){//all_mulike mass & high momentum
    pass_cut[10][0]=true;
  }

  if(kAllHist) MakeCutFlow();
  else MakeValidationPlot();

}

void OscNtupleManager::Process_pemumu(){

  //apply selection here

  if(!kDebugMode && wallv(0)<200) return;

  pass_cut[0][0]=true;
  if ( evis(0) > 30.0 && wall(0) > 200.0 && nhitac(0) <  nhitac_cut[skgen] ){//FC & FV
    pass_cut[1][0]=true;
  }

  //# of cherenkov ring
  if(nring(0)==2) pass_cut[2][0]=true;
  if(nring(0)==3) pass_cut[2][1]=true;
  if(nring(0)==4) pass_cut[2][2]=true;

  //#of mu-like ring
  if(n_mulike_angle==0) pass_cut[3][0]=true;
  if(n_mulike_angle==1) pass_cut[3][1]=true;
  if(n_mulike_angle==2) pass_cut[3][2]=true;
  if(n_mulike_angle==3) pass_cut[3][3]=true;

  // michel (decay) electron cut
  if(nDecayE==0) pass_cut[4][0]=true;
  if(nDecayE==1) pass_cut[4][1]=true;
  if(nDecayE==2) pass_cut[4][2]=true;
  if(nDecayE==3) pass_cut[4][3]=true;


  vector<TLorentzVector> mu_cand;
  TLorentzVector e_cand;
  mu_cand.clear();
  for(int r=0;r<nRing;r++){
    if( ( sqrt(fabs(probms(r,1)))-sqrt(fabs(probms(r,2))) ) < 0) e_cand = GetTLorentzVectorRing(r,1);
    else mu_cand.push_back(GetTLorentzVectorRing(r,2));
  }

  pass_cut[5][0]=true;//no pi0 mass cut


  TLorentzVector total_vec;
  if(n_mulike_angle==2) total_vec = mu_cand[0] + mu_cand[1] + e_cand;
  if(n_mulike_angle==1) total_vec = mu_cand[0] + e_cand;
  total_mass = total_vec.M();
  total_mom = total_vec.P();
  if(total_mass>800 && total_mass<1050 && total_mom<100){//total mass & low momentum
    pass_cut[6][0]=true;
  }
  if(total_mass>800 && total_mass<1050 && total_mom>100 && total_mom<250){//total mass & high momentum
    pass_cut[7][0]=true;
  }

  if(nNeutron==0) pass_cut[8][0]=true;

  TLorentzVector all_ring_vec;
  for(int r=0;r<nRing;r++){
    TLorentzVector temp_vec = 
    ( ( sqrt(fabs(probms(r,1)))-sqrt(fabs(probms(r,2))) ) < 0) ? GetTLorentzVectorRing(r,0) : GetTLorentzVectorRing(r,2);
    all_ring_vec += temp_vec;
  }
  all_ring_mass = all_ring_vec.M();
  all_ring_mom = all_ring_vec.P();

  if(kAllHist) MakeCutFlow();
  else MakeValidationPlot();

}

void OscNtupleManager::Process_pmuee(){

  //apply selection here

  if(!kDebugMode && wallv(0)<200) return;

  pass_cut[0][0]=true;
  if ( evis(0) > 30.0 && wall(0) > 200.0 && nhitac(0) <  nhitac_cut[skgen] ){//FC & FV
    pass_cut[1][0]=true;
  }

  //# of cherenkov ring
  if(nring(0)==2) pass_cut[2][0]=true;
  if(nring(0)==3) pass_cut[2][1]=true;
  //if(nring(0)==4) pass_cut[2][2]=true;
  if(nring(0)==2 || nring(0)==3) pass_cut[2][2]=true;

  //#of mu-like ring
  if(n_mulike_angle==0) pass_cut[3][0]=true;
  if(n_mulike_angle==1) pass_cut[3][1]=true;
  if(n_mulike_angle==2) pass_cut[3][2]=true;
  if(n_mulike_angle==3) pass_cut[3][3]=true;

  // michel (decay) electron cut
  if(nDecayE==0) pass_cut[4][0]=true;
  if(nDecayE==1) pass_cut[4][1]=true;
  if(nDecayE==2) pass_cut[4][2]=true;
  if(nDecayE==3) pass_cut[4][3]=true;

  vector<TLorentzVector> e_cand;
  TLorentzVector mu_cand;
  e_cand.clear();
  for(int r=0;r<nRing;r++){
    if( ( sqrt(fabs(probms(r,1)))-sqrt(fabs(probms(r,2))) ) < 0) e_cand.push_back(GetTLorentzVectorRing(r,1));
    else mu_cand = GetTLorentzVectorRing(r,2);
  }

  pass_cut[5][0]=true;//no pi0 mass cut


  TLorentzVector total_vec;
  if(n_elike_angle==2) total_vec = e_cand[0] + e_cand[1] + mu_cand;
  if(n_elike_angle==1) total_vec = e_cand[0] + mu_cand;
  total_mass = total_vec.M();
  total_mom = total_vec.P();
  if(total_mass>800 && total_mass<1050 && total_mom<100){//total mass & low momentum
    pass_cut[6][0]=true;
  }
  if(total_mass>800 && total_mass<1050 && total_mom>100 && total_mom<250){//total mass & high momentum
    pass_cut[7][0]=true;
  }

  if(nNeutron==0) pass_cut[8][0]=true;

  TLorentzVector all_ring_vec;
  for(int r=0;r<nRing;r++){
    TLorentzVector temp_vec = 
    ( ( sqrt(fabs(probms(r,1)))-sqrt(fabs(probms(r,2))) ) < 0) ? GetTLorentzVectorRing(r,0) : GetTLorentzVectorRing(r,2);
    all_ring_vec += temp_vec;
  }
  all_ring_mass = all_ring_vec.M();
  all_ring_mom = all_ring_vec.P();

  if(kAllHist) MakeCutFlow();
  else MakeValidationPlot();

}

void OscNtupleManager::Process_peemu(){

  //apply selection here

  if(!kDebugMode && wallv(0)<200) return;

  pass_cut[0][0]=true;
  if ( evis(0) > 30.0 && wall(0) > 200.0 && nhitac(0) <  nhitac_cut[skgen] ){//FC & FV
    pass_cut[1][0]=true;
  }

  //# of cherenkov ring
  if(nring(0)==2) pass_cut[2][0]=true;
  if(nring(0)==3) pass_cut[2][1]=true;
  //if(nring(0)==4) pass_cut[2][2]=true;
  if(nring(0)==2 || nring(0)==3) pass_cut[2][2]=true;

  //#of mu-like ring
  if(n_mulike_angle==0) pass_cut[3][0]=true;
  if(n_mulike_angle==1) pass_cut[3][1]=true;
  if(n_mulike_angle==2) pass_cut[3][2]=true;
  if(n_mulike_angle==3) pass_cut[3][3]=true;

  // michel (decay) electron cut
  if(nDecayE==0) pass_cut[4][0]=true;
  if(nDecayE==1) pass_cut[4][1]=true;
  if(nDecayE==2) pass_cut[4][2]=true;
  if(nDecayE==3) pass_cut[4][3]=true;

  vector<TLorentzVector> e_cand;
  TLorentzVector mu_cand;
  e_cand.clear();
  for(int r=0;r<nRing;r++){
    if( ( sqrt(fabs(probms(r,1)))-sqrt(fabs(probms(r,2))) ) < 0) e_cand.push_back(GetTLorentzVectorRing(r,1));
    else mu_cand = GetTLorentzVectorRing(r,2);
  }

  pass_cut[5][0]=true;//no pi0 mass cut


  TLorentzVector total_vec;
  if(n_elike_angle==2) total_vec = e_cand[0] + e_cand[1] + mu_cand;
  if(n_elike_angle==1) total_vec = e_cand[0] + mu_cand;
  total_mass = total_vec.M();
  total_mom = total_vec.P();
  if(total_mass>800 && total_mass<1050 && total_mom<100){//total mass & low momentum
    pass_cut[6][0]=true;
  }
  if(total_mass>800 && total_mass<1050 && total_mom>100 && total_mom<250){//total mass & high momentum
    pass_cut[7][0]=true;
  }

  if(nNeutron==0) pass_cut[8][0]=true;

  TLorentzVector all_ring_vec;
  for(int r=0;r<nRing;r++){
    TLorentzVector temp_vec = 
    ( ( sqrt(fabs(probms(r,1)))-sqrt(fabs(probms(r,2))) ) < 0) ? GetTLorentzVectorRing(r,0) : GetTLorentzVectorRing(r,2);
    all_ring_vec += temp_vec;
  }
  all_ring_mass = all_ring_vec.M();
  all_ring_mom = all_ring_vec.P();

  if(kAllHist) MakeCutFlow();
  else MakeValidationPlot();

}

void OscNtupleManager::ZeroStructure()
{

  os->Zero();
  event_type=-1;
  true_mode=-1;
  nRing = nring(0);
  nPar = npar(0);
  nPar2 = npar2(0);
  is_free_proton = (pmomv(0)>0.01)? 0 : 1;
  if(process_input=="fcmc" || process_input=="fcdt") is_free_proton = 0;
  closest_mass_pi0_reco=0.;
  total_mass=0.;
  total_mom=0.;
  all_ring_mass=0.;
  all_ring_mom=0.;

  //cout e-like or mu-like ring
  n_elike=0;n_elike_pattern=0;n_elike_angle=0;
  n_mulike=0;n_mulike_pattern=0;n_mulike_angle=0;
  for(int r=0;r<nRing;r++){
    float prob_angle = sqrt(fabs(probms(r,1)))-sqrt(fabs(probms(r,2)));
    float prob_pattern = prmslg(r,1) - prmslg(r,2);
    if(kDebugMode) {
      cout << "prob angle/pattern=" << prob_angle << "/" << prob_pattern << endl;
    }
    if(ip(r)==2) n_elike++;
    else if (ip(r)==3) n_mulike++;
    if(prob_angle<0) n_elike_angle++;
    else n_mulike_angle++;
    if(prob_pattern<0) n_elike_pattern++;
    else n_mulike_pattern++;
  }

  //count decay electron
  float llelectron = llm->GetMGMRELikelihood();
  nDecayE    = llm->GetDecayE();  // must be called after llBuild

  //neutron tagging
  if(skgen == SK4) nNeutron = ntag_nn(0);//only for sk4
  else nNeutron = 0;

  for(int c=0;c<10;c++) 
    for(int c2=0;c2<4;c2++)
      pass_cut[c][c2]=false;

  return;
}

bool OscNtupleManager::SetEventType( std::string mode , int &type)
{

  dataflag = 0; // default to MC processing
  std::string::size_type loc = mode.find("dt", 0 );
  if( loc != std::string::npos ) dataflag = 1; 


  if( mode == "fcdt" || mode == "fcmc" )
    return SetEventTypeFC( type );

  else 
  { 
    std::cout << " Error - unsupported mode : " << mode  << " exiting " << std::endl;
    exit(-1);
  }

  return false; //program flow should never make it here

}



bool OscNtupleManager::SetEventTypeFC( int &type)
{

  llMGMRE       = -1e7 ;
  llMGMRENUE    = -1e7 ;
  llMGMRENUEBar = -1e7 ;
  pi0ll         = -1e7 ;

  type = -1;
  bool fc_flag = false; 

  // get rid of NC tau events
  if ( abs(mode(0)) > 30 && abs(ipnu(0)) == 16 ) return false; 


  if ( evis(0) > 30.0 && wall(0) > 200.0 && nhitac(0) <  nhitac_cut[skgen] )
    fc_flag = true;

  // this is not a fc event
  if( ! fc_flag )
  {
    return false; 
  }

  float llelectron = llm->GetMGMRELikelihood();
  nDecayE    = llm->GetDecayE();  // must be called after llBuild

  llmpid  = llm->GetLLmpid() ; 
  llmue   = llm->GetLLmue () ; 
  llfmom  = llm->GetLLfmom() ; 
  lldpos  = llm->GetLLdpos() ; 


  float llnue    = llm->GetMGMRELikelihood_nue();  
  float llnuebar = llm->GetMGMRELikelihood_nuebar();


  bool multi_gev_flag = false;
  if( evis(0) >= 1330.0 ) multi_gev_flag = true;

  if( skgen == SK4 ){
    //os->nn = ntag_nn(0);//temporal
    //os->nn_mctruth=ntag_mctruth_nn(0);//temporal
  }

  ///////////////////
  // Sub-GeV Single Ring Selection 
  //

  // mu-like
  if( ! multi_gev_flag && nring(0) == 1 )
    if(  ip(0) == 3 && amomm(0) > 200.0 )
    {

      if( nDecayE == 0 ) 
      {

        type = SubGeV_mulike_0dcy; 

      }
      if( nDecayE == 1 )
      { 

        type = SubGeV_mulike_1dcy; 
      }
      if( nDecayE >= 2 )
      {
        type = SubGeV_mulike_2dcy; 
      }
      return true;
    }

  // e-like
  if( ! multi_gev_flag && nring(0) == 1 )
    if(  ip(0) == 2 && amome(0) > 100.0)
    { 
      pi0ll = llm->GetPi0Likelihood();

      //nDecayE was added by maggie in May 11
      if( (llm->GetPi0Selection() == 1 && nDecayE == 0) && ! ( kUseFiTQun &&  skgen == SK4 ) )
      {	     
        type = SubGeV_SingleRing_pi0like; 
        return true;
      }
      else
      {      

        // rvw hax -- skip fc taus! 20141114
        // needs to be removed after fiTQun is applied to tau MC
        if ( kUseFiTQun && skgen == SK4 && abs(ipnu(0)) != 16 ) 
        {

          if ( SubGeVElikeFQNCRejection() && nDecayE == 0 ) 
          {
            type = SubGeV_SingleRing_pi0like ; 
            return true  ; 
          }
        } 


        if( nDecayE == 0 )
        {	
          type = SubGeV_elike_0dcy; 
        }
        if( nDecayE >= 1 ) 
        {			
          type = SubGeV_elike_1dcy; 
        }
        return true;
      }
    }

  ///////////////////
  // Multi-GeV Single Ring Selection 
  //


  //multi GeV single-ring mu-like 
  if ( ip(0) == 3 && nring(0) == 1 && multi_gev_flag ) 
  {
    type = MultiGeV_mulike;
    return true;
  }       



  //multi GeV single-ring e-like 
  if ( ip(0) != 3 && nring(0) == 1 && multi_gev_flag ) 
  {    

    if(nDecayE > 0)
    {
      type = MultiGeV_elike_nue;		 
      return true;
    } 
    else if (nDecayE ==0)
    {     
      type = MultiGeV_elike_nuebar;
      return true;
    }

  }       
  ////////////////


  ///////////////////
  // Multi-GeV Multi Ring Selection 
  //

  // M.ost E.energetic R.ing index 
  int   mer=-1;
  float p;

  for( int j = 0 ; j < 3 ; j++ )
    ptot[j] = 0.;

  float pmax  = 0.; 
  for( int i = 0 ; i < nring(0) ; i++ )
  {
    p = ( prmslg(i,1) < prmslg(i,2) ? amome(i) : amomm(i) ); for( int j = 0 ; j < 3 ; j++ )
      ptot[j] += dir(i,j) * p;

    if( p > pmax ){ pmax = p; mer = i ;}  
  }

  // What is the identity of the 
  // M.ost E.nergetic R.ing
  int mer_id = ( prmslg(mer,1) < prmslg(mer,2) ? 1 : 2 ); // 1:e 2:mu


  // multi-ring mu-like sample
  if( mer_id  == 2 && nring(0) > 1) //mulike
    if( ( ! multi_gev_flag && evis(0) > 600.0 && pmax > 600.0 )  || multi_gev_flag )
    {
      type = MultiRing_mulike;
      return true;
    } 


  // MER is mulike events have already been selected so mer_id == elike cut
  // is reduntant here, but keep it for clarity
  if( multi_gev_flag && nring(0) > 1 && mer_id == 1  )
  {       

    llMGMRE = llelectron ;
    //std::cout << " MER ID: " << mer_id << " " << llelectron << std::endl;

    if( llelectron > 0 ) //e-like
    {
      llMGMRENUE    = llnue      ;
      llMGMRENUEBar = llnuebar   ;

      if( llnue > llnuebar )
      {         
        type= MultiRing_elike_nue;
        return true;
      }
      else
      {     
        type = MultiRing_elike_nuebar;  
        return true;
      }
    }
    else // things that failed the likelihood selection 
    {
      type = MultiRingOther_1;
      return true;
    }



  }// end multi-ting multi-gev chech

  ///////////////////
  // 2-ring Pi0-like 
  //
  float mass = sqrt(2.*amome(0)*amome(1)
      *(1.-msdir(0,0,1)*msdir(1,0,1)
        -msdir(0,1,1)*msdir(1,1,1)
        -msdir(0,2,1)*msdir(1,2,1) ));

  if( ! multi_gev_flag && nring(0) == 2 )
    if( ip(0) == 2 && ip(1) == 2 && nDecayE == 0 )
      if( mass > 85.0 && mass < 215.0) 
      {  
        type = SubGeV_pi0like; 
        return true;
      }


  // failed FC classification  

  return false;
}


bool OscNtupleManager::SubGeVElikeFQNCRejection()
{

  bool RejectAsNC = false ;
  if (     130.0 < ( fq1rnll(0,1) - fqpi0nll(0) )
      &&  95.0 < fqpi0mass(0) 
      && 185.0 > fqpi0mass(0) 
     )
    RejectAsNC = true ;

  return RejectAsNC;
}

float OscNtupleManager::CalcDrRingVector(int id_ring, int id_vector){
  float dx = dir(id_ring,0) - dirv(id_vector,0);
  float dy = dir(id_ring,1) - dirv(id_vector,1);
  float dz = dir(id_ring,2) - dirv(id_vector,2);
  float dr = sqrt(dx*dx+dy*dy+dz*dz);
  return dr;
}

float OscNtupleManager::CalcDrRingVector2(int id_ring, int id_vector2){
  float dx = dir(id_ring,0) - dirv2(id_vector2,0);
  float dy = dir(id_ring,1) - dirv2(id_vector2,1);
  float dz = dir(id_ring,2) - dirv2(id_vector2,2);
  float dr = sqrt(dx*dx+dy*dy+dz*dz);
  return dr;
}

float OscNtupleManager::CalcDrVectorVector(int id1_vector, int id2_vector){
  float dx = dirv(id1_vector,0) - dirv(id2_vector,0);
  float dy = dirv(id1_vector,1) - dirv(id2_vector,1);
  float dz = dirv(id1_vector,2) - dirv(id2_vector,2);
  float dr = sqrt(dx*dx+dy*dy+dz*dz);
  return dr;
}

float OscNtupleManager::CalcDrVector2Vector2(int id1_vector2, int id2_vector2){
  float dx = dirv2(id1_vector2,0) - dirv2(id2_vector2,0);
  float dy = dirv2(id1_vector2,1) - dirv2(id2_vector2,1);
  float dz = dirv2(id1_vector2,2) - dirv2(id2_vector2,2);
  float dr = sqrt(dx*dx+dy*dy+dz*dz);
  return dr;
}

TLorentzVector OscNtupleManager::GetTLorentzVectorVector(int index){
  TLorentzVector ans_vector;
  float mom = pmomv(index);
  float px = mom * dirv(index,0);//MeV
  float py = mom * dirv(index,1);//MeV
  float pz = mom * dirv(index,2);//MeV
  float mass = 0.;//MeV
  int pidgen = ipv(index);
  if(pidgen==1 || pidgen==4) mass = 0.;//gamma neutrino
  if(pidgen==2 || pidgen==3) mass = 0.511;//e+ e-
  if(pidgen==5 || pidgen==6) mass = 105.7;//mu+ mu-
  if(pidgen==7) mass = 135.;//pi0
  if(pidgen==8 || pidgen==9) mass = 139.6;//pi+ pi-
  if(pidgen==13) mass = 939.6;//neutron
  if(pidgen==14) mass = 938.3;//proton
  float energy = sqrt(mom*mom + mass*mass);
  ans_vector.SetPxPyPzE(px,py,pz,energy);
  return ans_vector;
}

TLorentzVector OscNtupleManager::GetTLorentzVectorVector2(int index){
  TLorentzVector ans_vector;
  float mom = pmomv2(index);
  float px = mom * dirv2(index,0);//MeV
  float py = mom * dirv2(index,1);//MeV
  float pz = mom * dirv2(index,2);//MeV
  float mass = 0.;//MeV
  int pidgen = ipv2(index);
  if(pidgen==1 || pidgen==4) mass = 0.;//gamma neutrino
  if(pidgen==2 || pidgen==3) mass = 0.511;//e+ e-
  if(pidgen==5 || pidgen==6) mass = 105.7;//mu+ mu-
  if(pidgen==7) mass = 135.;//pi0
  if(pidgen==8 || pidgen==9) mass = 139.6;//pi+ pi-
  if(pidgen==13) mass = 939.6;//neutron
  if(pidgen==14) mass = 938.3;//proton
  float energy = sqrt(mom*mom + mass*mass);
  ans_vector.SetPxPyPzE(px,py,pz,energy);
  return ans_vector;
}

TLorentzVector OscNtupleManager::GetTLorentzVectorRing(int index, int type){
  //type 0:gamma 1:electron 2:muon
  TLorentzVector ans_vector;
  float mome = amome(index);
  float momm = amomm(index);
  float px=0.,py=0.,pz=0.,mass=0.,energy=0.;
  if(type==0 || type==1){
    px = mome * msdir(index,0,1);//MeV
    py = mome * msdir(index,1,1);//MeV
    pz = mome * msdir(index,2,1);//MeV
    if(type==0) mass = 0.;//MeV
    if(type==1) mass = 0.511;//MeV
    energy = sqrt(mome*mome + mass*mass);
  }
  if(type==2){
    px = momm * msdir(index,0,2);//MeV
    py = momm * msdir(index,1,2);//MeV
    pz = momm * msdir(index,2,2);//MeV
    mass = 105.7;
    energy = sqrt(momm*momm + mass*mass);
  }
  ans_vector.SetPxPyPzE(px,py,pz,energy);
  return ans_vector;
}

TVector3 OscNtupleManager::GetTVectorVector(int index){
  TVector3 ans_vector;
  float x = dirv(index,0);
  float y = dirv(index,1);
  float z = dirv(index,2);
  ans_vector.SetXYZ(x,y,z);
  return ans_vector;
}

TVector3 OscNtupleManager::GetTVectorVector2(int index){
  TVector3 ans_vector;
  float x = dirv2(index,0);
  float y = dirv2(index,1);
  float z = dirv2(index,2);
  ans_vector.SetXYZ(x,y,z);
  return ans_vector;
}

TVector3 OscNtupleManager::GetTVectorRing(int index){
  TVector3 ans_vector;
  float x = dir(index,0);
  float y = dir(index,1);
  float z = dir(index,2);
  ans_vector.SetXYZ(x,y,z);
  return ans_vector;
}

TVector3 OscNtupleManager::GetTVectorRingPID(int index, int id){//1 e-like, 2 mu-like
  TVector3 ans_vector;
  float x = msdir(index,0,id);
  float y = msdir(index,1,id);
  float z = msdir(index,2,id);
  ans_vector.SetXYZ(x,y,z);
  return ans_vector;
}


double OscNtupleManager::Oscillate(){

  return Neighbor3D();

}

double OscNtupleManager::Neighbor3D(){

  // Converts PDG->internal OscType    e:12->1 , mu: 14->2, tau:16->3 //
  int    NuOscillatedTo = int( abs( ipnu(0) )/2 - 5 );
  //std::cout << "ipnu/NuOscillatedTo=" << ipnu(0) << "/" << NuOscillatedTo << std::endl;
  double Oscillated = 0;	
  double MCWeight = GetMCweight();
  double FactorE = 0, FactorMu = 0;
  int    NuStart;
  int    i,j;
  double Energy;
  bool   kUseSquaredThetas = true ;

  // Mass Hierarchy Flag 
  // ( correction for solar term is done internally )
  bool kInverted = false; //temporal
  double hFactor           = ( kInverted ? -1.0 : 1.0 ) ; 

  // for the fortran computation
  float  XXPATH[20];
  int    IDNU2    = NuOscillatedTo;
  float  DirNeu   = -1. * dirnu(0,2); 
  float  fortranE = pnu(0);
  double Path;

  int this_itype;
  bool sc = SetEventType(process_input, this_itype);
  this_itype = this_itype + (skgen * (EndOfTypes-1)) ;  
  //std::cout << "this_itype=" << this_itype << std::endl;
  //std::cout << "skgen/EndOfTypes=" << skgen << "/" << EndOfTypes << std::endl;

  BuildEnergyFriend(this_itype);
  // averaging counters
  double PathLengthAverage = 0.0;
  double FullPathAve = 0.0;
  double rms = ErmsHax(0) ;
  //double rms = 0;//temporal


  // Zero out any erroneou NC Tau events that 
  // may be included
  if ( abs( mode(0) ) >= 30 && abs(NuOscillatedTo)==3 )
    return 0;

  // NC events only recieve their solar flux weighting
  if ( abs( mode(0) ) >= 30 ) 
  {
    return MCWeight;
  }

  if ( ipnu(0) ==  12 ) IDNU2 = 1;
  if ( ipnu(0) == -12 ) IDNU2 = 3;
  if ( ipnu(0) ==  14 ) IDNU2 = 2;
  if ( ipnu(0) == -14 ) IDNU2 = 4;
  if ( ipnu(0) ==  16 ) IDNU2 = 2;
  if ( ipnu(0) == -16 ) IDNU2 = 4;

  Energy   = pnu(0);
  if( rms / Energy  >= 0.5 ) rms = 0.5 * Energy ;

  // Flux factor adjustment : 
  // if nue      ,  FactorE   = 1 , otherwise reweight from muonflux to electron flux
  // if numu/tau ,  FactorMu  = 1 , otherwise reweight from electron to muon flux 
  FactorE  = ( NuOscillatedTo == 1 ?  1 : GetHondaFluxRatio( 2 ) );  
  FactorMu = ( NuOscillatedTo == 2 || NuOscillatedTo == 3 ? 1 : GetHondaFluxRatio( 1 ) );
  //std::cout << "FactorE/FactorMu=" << FactorE << "/" << FactorMu << std::endl;

  int NEAve = nEAveHax(0);
  //int NEAve = 0;//temporal


  double EAverages [5] =  { Energy           ,
    Energy - rms     , 
    Energy - rms *0.5 , 
    Energy + rms *0.5 , 
    Energy + rms      
  }; 


  if( NEAve == 0 )
  {
    EAverages[1] = 1.50 * Energy ;
    EAverages[2] = 1.25 * Energy ;
    EAverages[3] = 0.75 * Energy ;
    EAverages[4] = 0.50 * Energy ;
  }

  float SKheight = 0.380;
  m_fort->nebaseline(XXPATH, IDNU2, DirNeu, fortranE, SKheight);

  for( j = 0 ; j < 20 ; j++)
    FullPathAve += (double) XXPATH[j] /  20.0;     

  //std::cout << "FullPathAve=" << FullPathAve << std::endl;

  double AvePath   [5] ;
  AvePath [0]  =   FullPathAve * 1.0e5;  // in [cm]

  // These are lifted from osc3d_map.F (arbitrary assignment)
  int PathAveIndices[4][4] = {  // Original Fortran Indexing 
    10 , 9 , 11 , 12  , 
    8 , 7 , 13 , 14  ,
    6 , 5 , 15 , 16  ,
    4 , 3 , 17 , 18 
  };

  for( i = 0 ; i < 4 ; i++ ) 
  {
    AvePath [i+1]  = 0.0 ;
    for( j = 0 ; j < 4 ; j++ ) 
    {
      // convert to c-indexing
      int index = PathAveIndices[i][j] - 1;
      AvePath [i+1]  +=   (double) XXPATH[ index ] / 4.0   ; 
    }
    AvePath [i+1] *= 1.0e5  ; // convert to cm
  }

  // Set once with arbitrary production height,
  // will be overridden with call to SetPathLength below
  // false argument stops matter profile specification 
  bNu->DefinePath( dirnu(0,2), 15.00, false ); 

  //bNu->SetMNS( P.Get("S12") , P.Get("S13") , P.Get("S23") , 
  //             P.Get("M12") ,      hFactor * P.Get("M23") , 
  //             P.Get("CP")  , 
  //             Energy       , kUseSquaredThetas ,  E->GetPDG() );  

  bNu->SetMNS( 0.309 , 0.0219 , 0.3 , 
      0.0000765 ,      hFactor * 0.0025 , 
      4.19  , 
      Energy       , kUseSquaredThetas ,  ipnu(0) );  


  int lnPaths = 5; 
  int Type = (this_itype - 1) % global::nEventTypes;
  //std::cout << "Type=" << Type << std::endl;

  // skip energy averaging for 
  // upthrough types
  if( Type == global::UpThruNonShower_mu ||
      Type == global::UpThruShower_mu  )  
    lnPaths = 1 ;


  double PathLength;
  double probe = 0.;
  double probm = 0.;

  for (i = 0; i < lnPaths ; i++ ) // energy path length aver
  {
    PathLengthAverage += 1.00;                                                                                                                       
    bNu->SetPathLength( AvePath[i]    );	 
    bNu->SetEnergy    ( EAverages[i]  );

    NuStart = 1;
    NuStart *= ( ipnu(0) < 0 ? -1 : 1 );

    // Only need to propage once 
    bNu->propagate( NuStart );

    Oscillated += bNu->GetProb( NuStart, NuOscillatedTo )*MCWeight*FactorE;
    probe       = bNu->GetProb( NuStart, NuOscillatedTo );

    NuStart = 2;
    NuStart *= ( ipnu(0) < 0 ? -1 : 1 );

    Oscillated += bNu->GetProb( NuStart, NuOscillatedTo )*MCWeight*FactorMu;
    probm       = bNu->GetProb( NuStart, NuOscillatedTo );

  } // End of Pathlength averaging

  Oscillated /= ( PathLengthAverage );

  if(isnan(Oscillated)) Oscillated = 0.;

  return Oscillated;

}

double OscNtupleManager::GetMCweight(){
  //std::cout << "OscNtupleManager::GetMCweight" << std::endl;

  double weightx=0.;

  if( skgen == SK1 )
  {
    //sk1
    weightx = (0.7*flxh11(0)+0.3*flxh11(2))/flxh06(1);    
  }

  if( skgen == SK2 )
  {
    weightx = (0.3*flxh11(0)+0.7*flxh11(2))/flxh06(1);
  }

  // NB this will have to be updated for 
  // for SK-IV
  if( skgen == SK3 )
  {
    weightx = (1.0*flxh11(0)+0.0*flxh11(2))/flxh06(1);
  }  

  // NB this will have to be updated for 
  // for SK-IV
  if( skgen >= SK4 )
  {
    weightx = (0.65*flxh11(0)+0.35*flxh11(2))/flxh06(1);
  }  

  if ( abs(ipnu(0)) == 16 ) weightx = 1.0; 
  //std::cout << "skgen/weightx=" << skgen << "/" << weightx << std::endl;

  return weightx;

}

double OscNtupleManager::GetHondaFluxRatio( int NuType )
{
  //std::cout << "GetHondaFluxRatio" << std::endl;
  float Solar    = 0.5;
  float d_pnu = pnu(0);
  float d_dirnu[3];
  for( int i = 0 ; i < 3 ; i++ ){
    d_dirnu[i] = dirnu(0,i);
    //std::cout << "d_dirnu=" << d_dirnu[i] << std::endl;
  }

  int   type[] = {12, 14};
  double flxho[2];
  for( int i = 0 ; i < 2 ; i++ )
  {
    int nu = type[i] * ( ipnu(0) < 1 ? -1 : 1 );
    //std::cout << "d_pnu/d_dir_nu/Solar/nu=" << d_pnu << "/" << d_dirnu << "/" << Solar << "/" << nu << std::endl;
    //flxho[i] = fnhonfx11_( &d_pnu, d_dirnu, &Solar, &nu ) ;
    flxho[i] = m_fort->calc_flux( d_pnu, d_dirnu, Solar, nu ) ;
    //std::cout << "flxho=" << flxho[i] << std::endl;
  }

  if ( NuType == 1 )   
    return (double) flxho[1] / flxho[0];

  if ( NuType == 2 || NuType == 3 )
    return (double) flxho[1] / flxho[0];

  std::cerr << "Returning 0 from SKEventParser::GetHondaFluxRatio for type " << NuType << std::endl;
  return 0;

}

void OscNtupleManager::BuildEnergyFriend(int type){
  //std::cout << "OscNtupleManager::BuildEnergyFriend" << std::endl;
  int zbin = GetZBin(-1.*dir(2,0),type);
  //std::cout << "zbin=" << zbin << std::endl;


}

int OscNtupleManager::GetZBin(float zenith, int type){
  float ncz = 10.; 
  float czstart = -1.;
  float czstop =1.; 

  if( zenith >= 1.0 ) zenith = 0.9999999;

  int cmp  =  ( type -1 ) % 19 ;
  cmp += 1; 

  // print itype, OscTypes.nTypes,  cmp , OscTypes.SubGeV_SingleRing_pi0like 
  if ( cmp == SubGeV_SingleRing_pi0like  ) ncz = 1.0 ;
  if ( cmp == SubGeV_elike_1dcy          ) ncz = 1.0 ;
  if ( cmp == SubGeV_mulike_2dcy         ) ncz = 1.0 ; 
  if ( cmp == SubGeV_pi0like             ) ncz = 1.0 ; 

  if ( cmp == UpStop_mu || cmp == UpThruNonShower_mu || cmp == UpThruShower_mu     ) czstop = 0.0 ; 

  float binWidth = (czstop - czstart)/ncz ;
  int bin      = int( (zenith - czstart) /binWidth ); 

  return bin;
}

void OscNtupleManager::FillNtuple(){

  int   type;
  SetEventType(process_input,type);

  // convert to indexing for other 
  // skgenerations 
  o_itype = type + (skgen * (EndOfTypes-1)) ;  
  o_ipnu = ipnu(0);
  o_pnu  = pnu(0);
  o_mode = mode(0);


  //  Single-Ring Fully Contained Events
  if( type <  MultiRing_elike_nue  ){
    for( int i=0 ; i < 3 ; i++ ) o_dir[i] = dir(0,i);
    // are we an e-like type?
    if( probms(0,1) - probms(0,2) > 0 ) o_amom = amome(0);
    else o_amom = amomm(0);
  }

  //  Multi-Ring Fully Contained Events
  if( type == MultiRing_elike_nue || type == MultiRing_elike_nuebar || type == MultiRingOther_1 ) {
    for( int i=0 ; i < 3 ; i++ ) o_dir[i] = dirtotmue(i);
    o_amom = etotmue(0);
  }  

  if (type == MultiRing_mulike) {
    // ptot is computed below in 
    // SetEventTypeFC
    float norm  = ptot[0]*ptot[0];
    norm += ptot[1]*ptot[1];
    norm += ptot[2]*ptot[2];
    for( int i=0 ; i < 3 ; i++ )
      o_dir[i] = ptot[i]/sqrt(norm);
    o_amom = evis(0);
  }

  otree->Fill();

}

void OscNtupleManager::MakeOscillationPlot(){
  if(kDebugMode){
    std::cout << "interacion_type=" << interaction_type << std::endl;
    std::cout << "weight mc/osc=" << mc_weight << "/" << osc_weight << std::endl;
  }
  if(interaction_type!=-1)  {
    m_hSvc.h1D(Form("nRing_type%d",interaction_type),"","")->Fill(nring(0),weight);
    m_hSvc.h1D(Form("nRing_type%d_osc",interaction_type),"","")->Fill(nring(0),weight*osc_weight);
  }
  if(interaction_type>=1 && interaction_type<=10 && interaction_type!=7){//single ring FC event
    //std::cout << "single ring event?" << std::endl;
    //std::cout << "nring=" << nring(0) << std::endl;
    float zenith_angle = -1.*dir(0,2);
    float momentum = (ip(0)==2)? amome(0) : amomm(0);
    m_hSvc.h1D(Form("zenith_angle_type%d",interaction_type),"","")->Fill(zenith_angle,weight);
    m_hSvc.h1D(Form("zenith_angle_type%d_osc",interaction_type),"","")->Fill(zenith_angle,weight*osc_weight);
    m_hSvc.h1D(Form("momentum_log10_type%d",interaction_type),"","")->Fill(log10(momentum),weight);
    m_hSvc.h1D(Form("momentum_log10_type%d_osc",interaction_type),"","")->Fill(log10(momentum),weight*osc_weight);
    for(int m=0;m<n_range_momentum-1;m++){
      if(momentum>range_momentum[m] && momentum<range_momentum[m+1]){
        m_hSvc.h1D(Form("zenith_angle_type%d_mom%d_%d",interaction_type,range_momentum[m],range_momentum[m+1]),"","")->Fill(zenith_angle,weight);
        m_hSvc.h1D(Form("zenith_angle_type%d_mom%d_%d_osc",interaction_type,range_momentum[m],range_momentum[m+1]),"","")->Fill(zenith_angle,weight*osc_weight);
      }
    }

  }
  else if ((interaction_type>=11 && interaction_type<=14) || interaction_type==7){//multi ring FC
    int id_energetic=-1;
    float max_mom=0;
    for(int r=0;r<nRing;r++){
      float momentum = (ip(r)==2)? amome(r) : amomm(r);
      if(momentum>max_mom){
        max_mom=momentum;
        id_energetic=r;
      }
    }
    for(int m=0;m<n_range_momentum-1;m++){
      if(max_mom>range_momentum[m] && max_mom<range_momentum[m+1]){
        float zenith_angle = -1.*dir(id_energetic,2);
        m_hSvc.h1D(Form("zenith_angle_type%d_mom%d_%d",interaction_type,range_momentum[m],range_momentum[m+1]),"","")->Fill(zenith_angle,weight);
        m_hSvc.h1D(Form("zenith_angle_type%d_mom%d_%d_osc",interaction_type,range_momentum[m],range_momentum[m+1]),"","")->Fill(zenith_angle,weight*osc_weight);
      }
    }
  }
}

void OscNtupleManager::MakeCutFlow(){

  for(int r=0;r<3;r++){// # of ring
    for(int mu=0;mu<2;mu++){// # of mu-like ring
      for(int m=0;m<2;m++){//# of michel electron
        for(int c=0;c<8;c++){//cut
          if(c<6){
            if(!pass_cut[c][0] && c!=2 && c!=3 && c!=4)  break;
            if((!pass_cut[c][r] && c==2) || (!pass_cut[c][mu] && c==3) || (!pass_cut[c][m] && c==4) ) break;
            MakeBasicPlot(c,r,mu,m);
          }
          else if(pass_cut[c][0]) {
            MakeBasicPlot(c,r,mu,m);
            if(pass_cut[8][0]) MakeBasicPlot(c+2,r,mu,m);//neutron tagging
          }
        }
      }
    }
  }

}

void OscNtupleManager::MakeCutFlowValidate(){

  for(int c=0;c<11;c++){//cut
    if(c<7 && !pass_cut[c][0])  break;
    else if(!pass_cut[c][0]) continue; 
    MakeBasicPlot(c,0,0,0);
  }

}

void OscNtupleManager::MakeBasicPlot(int c, int r, int mu, int p){//cut #, michel e cut, 


  if(process_input!="fcdt" || total_mass<800 || total_mass>1050 || total_mom>250){//for blind analysis
    m_hSvc.h1D(Form("cut_flow_nring%d_mulike%d_michel%d",r,mu,p),"","")->Fill(c,weight*osc_weight);
    m_hSvc.h1D(Form("cut_flow_nring%d_mulike%d_michel%d_fp%d",r,mu,p,is_free_proton),"","")->Fill(c,weight*osc_weight);
  }
  m_hSvc.h1D(Form("interaction_mode_cut%d_nring%d_mulike%d_michel%d",c,r,mu,p),"","")->Fill(true_mode,weight*osc_weight);
  m_hSvc.h1D(Form("distance_to_wall_cut%d_nring%d_mulike%d_michel%d",c,r,mu,p),"","")->Fill(wall(0),weight*osc_weight);
  m_hSvc.h1D(Form("visible_energy_cut%d_nring%d_mulike%d_michel%d",c,r,mu,p),"","")->Fill(evis(0),weight*osc_weight);
  m_hSvc.h1D(Form("nhit_OD_cut%d_nring%d_mulike%d_michel%d",c,r,mu,p),"","")->Fill(nhitac(0),weight*osc_weight);
  m_hSvc.h1D(Form("nRing_cut%d_nring%d_mulike%d_michel%d",c,r,mu,p),"","")->Fill(nRing,weight*osc_weight);
  m_hSvc.h1D(Form("nElikeRing_angle_nring%d_cut%d_nring%d_mulike%d_michel%d",nRing,c,r,mu,p),"","")->Fill(n_elike_angle,weight*osc_weight);
  m_hSvc.h1D(Form("nElikeRing_pattern_nring%d_cut%d_nring%d_mulike%d_michel%d",nRing,c,r,mu,p),"","")->Fill(n_elike_pattern,weight*osc_weight);
  m_hSvc.h1D(Form("nMulikeRing_nring%d_cut%d_nring%d_mulike%d_michel%d",nRing,c,r,mu,p),"","")->Fill(n_mulike,weight*osc_weight);
  m_hSvc.h1D(Form("nMulikeRing_angle_nring%d_cut%d_nring%d_mulike%d_michel%d",nRing,c,r,mu,p),"","")->Fill(n_mulike_angle,weight*osc_weight);
  m_hSvc.h1D(Form("nMulikeRing_pattern_nring%d_cut%d_nring%d_mulike%d_michel%d",nRing,c,r,mu,p),"","")->Fill(n_mulike_pattern,weight*osc_weight);
  m_hSvc.h1D(Form("n_michel_electron_cut%d_nring%d_mulike%d_michel%d",c,r,mu,p),"","")->Fill(nDecayE,weight*osc_weight);
  m_hSvc.h1D(Form("ntag_multiplicity_cut%d_nring%d_mulike%d_michel%d",c,r,mu,p),"","")->Fill(nNeutron,weight*osc_weight);
  m_hSvc.h1D(Form("mass_pi0_reco_elike%d_cut%d_nring%d_mulike%d_michel%d",n_elike_pattern,c,r,mu,p),"","")->Fill(closest_mass_pi0_reco,weight*osc_weight);

  m_hSvc.h1D(Form("interaction_mode_cut%d_nring%d_mulike%d_michel%d_fp%d",c,r,mu,p,is_free_proton),"","")->Fill(true_mode,weight*osc_weight);
  m_hSvc.h1D(Form("distance_to_wall_cut%d_nring%d_mulike%d_michel%d_fp%d",c,r,mu,p,is_free_proton),"","")->Fill(wall(0),weight*osc_weight);
  m_hSvc.h1D(Form("visible_energy_cut%d_nring%d_mulike%d_michel%d_fp%d",c,r,mu,p,is_free_proton),"","")->Fill(evis(0),weight*osc_weight);
  m_hSvc.h1D(Form("nhit_OD_cut%d_nring%d_mulike%d_michel%d_fp%d",c,r,mu,p,is_free_proton),"","")->Fill(nhitac(0),weight*osc_weight);
  m_hSvc.h1D(Form("nRing_cut%d_nring%d_mulike%d_michel%d_fp%d",c,r,mu,p,is_free_proton),"","")->Fill(nRing,weight*osc_weight);
  m_hSvc.h1D(Form("nElikeRing_angle_nring%d_cut%d_nring%d_mulike%d_michel%d_fp%d",nRing,c,r,mu,p,is_free_proton),"","")->Fill(n_elike_angle,weight*osc_weight);
  m_hSvc.h1D(Form("nElikeRing_pattern_nring%d_cut%d_nring%d_mulike%d_michel%d_fp%d",nRing,c,r,mu,p,is_free_proton),"","")->Fill(n_elike_pattern,weight*osc_weight);
  m_hSvc.h1D(Form("nMulikeRing_angle_nring%d_cut%d_nring%d_mulike%d_michel%d_fp%d",nRing,c,r,mu,p,is_free_proton),"","")->Fill(n_mulike_angle,weight*osc_weight);
  m_hSvc.h1D(Form("nMulikeRing_pattern_nring%d_cut%d_nring%d_mulike%d_michel%d_fp%d",nRing,c,r,mu,p,is_free_proton),"","")->Fill(n_mulike_pattern,weight*osc_weight);
  m_hSvc.h1D(Form("n_michel_electron_cut%d_nring%d_mulike%d_michel%d_fp%d",c,r,mu,p,is_free_proton),"","")->Fill(nDecayE,weight*osc_weight);
  m_hSvc.h1D(Form("ntag_multiplicity_cut%d_nring%d_mulike%d_michel%d_fp%d",c,r,mu,p,is_free_proton),"","")->Fill(nNeutron,weight*osc_weight);
  m_hSvc.h1D(Form("mass_pi0_reco_elike%d_cut%d_nring%d_mulike%d_michel%d_fp%d",n_elike_pattern,c,r,mu,p,is_free_proton),"","")->Fill(closest_mass_pi0_reco,weight*osc_weight);
  if(process_input!="fcdt" || total_mass<800 || total_mass>1050){//for blind analysis
    m_hSvc.h1D(Form("mass_proton_reco_cut%d_nring%d_mulike%d_michel%d",c,r,mu,p),"","")->Fill(total_mass,weight*osc_weight);
    m_hSvc.h1D(Form("mass_proton_reco_cut%d_nring%d_mulike%d_michel%d_fp%d",c,r,mu,p,is_free_proton),"","")->Fill(total_mass,weight*osc_weight);
    m_hSvc.h1D(Form("mass_all_ring_reco_nring%d_cut%d_nring%d_mulike%d_michel%d_fp%d",nRing,c,r,mu,p,is_free_proton),"","")->Fill(all_ring_mass,weight*osc_weight);
    m_hSvc.h1D(Form("mass_all_mulike_reco_nring%d_cut%d_nring%d_mulike%d_michel%d_fp%d",nRing,c,r,mu,p,is_free_proton),"","")->Fill(all_mulike_mass,weight*osc_weight);
  }
  if(process_input!="fcdt" || total_mom>250){//for blind analysis
    m_hSvc.h1D(Form("mom_proton_reco_cut%d_nring%d_mulike%d_michel%d",c,r,mu,p),"","")->Fill(total_mom,weight*osc_weight);
    m_hSvc.h1D(Form("mom_proton_reco_cut%d_nring%d_mulike%d_michel%d_fp%d",c,r,mu,p,is_free_proton),"","")->Fill(total_mom,weight*osc_weight);
    m_hSvc.h1D(Form("mom_all_ring_reco_nring%d_cut%d_nring%d_mulike%d_michel%d_fp%d",nRing,c,r,mu,p,is_free_proton),"","")->Fill(all_ring_mom,weight*osc_weight);
    m_hSvc.h1D(Form("mom_all_mulike_reco_nring%d_cut%d_nring%d_mulike%d_michel%d_fp%d",nRing,c,r,mu,p,is_free_proton),"","")->Fill(all_mulike_mom,weight*osc_weight);
  }
  if(process_input!="fcdt" || total_mass<800 || total_mass>1050 || total_mom>250){//for blind analysis
    m_hSvc.graph(Form("mass_mom_proton_reco_cut%d_nring%d_mulike%d_michel%d",c,r,mu,p),"","")->SetPoint(graph_point_2,total_mass,total_mom);
    m_hSvc.graph(Form("all_ring_mass_mom_proton_reco_cut%d_nring%d_mulike%d_michel%d",c,r,mu,p),"","")->SetPoint(graph_point_2,all_ring_mass,all_ring_mom);
    m_hSvc.graph(Form("all_mulike_mass_mom_proton_reco_cut%d_nring%d_mulike%d_michel%d",c,r,mu,p),"","")->SetPoint(graph_point_2,all_mulike_mass,all_mulike_mom);
    m_hSvc.graph(Form("mass_mom_proton_reco_cut%d_nring%d_mulike%d_michel%d_fp%d",c,r,mu,p,is_free_proton),"","")->SetPoint(graph_point_2,total_mass,total_mom);
    m_hSvc.graph(Form("all_ring_mass_mom_proton_reco_cut%d_nring%d_mulike%d_michel%d_fp%d",c,r,mu,p,is_free_proton),"","")->SetPoint(graph_point_2,all_ring_mass,all_ring_mom);
    m_hSvc.graph(Form("all_mulike_mass_mom_proton_reco_cut%d_nring%d_mulike%d_michel%d_fp%d",c,r,mu,p,is_free_proton),"","")->SetPoint(graph_point_2,all_mulike_mass,all_mulike_mom);
  }
  if(kCheckBkg){//decault is fcmc
    m_hSvc.h1D(Form("interaction_mode_cut%d_nring%d_mulike%d_michel%d_mode%d",c,r,mu,p,true_mode),"","")->Fill(true_mode,weight*osc_weight);
    m_hSvc.h1D(Form("distance_to_wall_cut%d_nring%d_mulike%d_michel%d_mode%d",c,r,mu,p,true_mode),"","")->Fill(wall(0),weight*osc_weight);
    m_hSvc.h1D(Form("visible_energy_cut%d_nring%d_mulike%d_michel%d_mode%d",c,r,mu,p,true_mode),"","")->Fill(evis(0),weight*osc_weight);
    m_hSvc.h1D(Form("nhit_OD_cut%d_nring%d_mulike%d_michel%d_mode%d",c,r,mu,p,true_mode),"","")->Fill(nhitac(0),weight*osc_weight);
    m_hSvc.h1D(Form("nRing_cut%d_nring%d_mulike%d_michel%d_mode%d",c,r,mu,p,true_mode),"","")->Fill(nRing,weight*osc_weight);
    m_hSvc.h1D(Form("nElikeRing_nring%d_cut%d_nring%d_mulike%d_michel%d_mode%d",nRing,c,r,mu,p,true_mode),"","")->Fill(n_elike,weight*osc_weight);
    m_hSvc.h1D(Form("n_michel_electron_cut%d_nring%d_mulike%d_michel%d_mode%d",c,r,mu,p,true_mode),"","")->Fill(nDecayE,weight*osc_weight);
    m_hSvc.h1D(Form("ntag_multiplicity_cut%d_nring%d_mulike%d_michel%d_mode%d",c,r,mu,p,true_mode),"","")->Fill(nNeutron,weight*osc_weight);
    m_hSvc.h1D(Form("mass_pi0_reco_elike%d_cut%d_nring%d_mulike%d_michel%d_mode%d",n_elike_pattern,c,r,mu,p,true_mode),"","")->Fill(closest_mass_pi0_reco,weight*osc_weight);
    m_hSvc.h1D(Form("mass_proton_reco_cut%d_nring%d_mulike%d_michel%d_mode%d",c,r,mu,p,true_mode),"","")->Fill(total_mass,weight*osc_weight);
    m_hSvc.h1D(Form("mom_proton_reco_cut%d_nring%d_mulike%d_michel%d_mode%d",c,r,mu,p,true_mode),"","")->Fill(total_mom,weight*osc_weight);
    m_hSvc.graph(Form("mass_mom_proton_reco_cut%d_nring%d_mulike%d_michel%d_mode%d",c,r,mu,p,true_mode),"","")->SetPoint(graph_point_2,total_mass,total_mom);
  }
  graph_point_2++;

}

void OscNtupleManager::MakeValidationPlot(){

  //if ( evis(0) > 30.0 && wall(0) > 200.0 && nhitac(0) <  nhitac_cut[skgen] ){}//FC & FC cut
  //else return;

  m_hSvc.graph("tgraph_test","","")->SetPoint(graph_point,graph_point,graph_point);

  float vertex_r_ring = sqrt(pos(0)*pos(0)+pos(1)*pos(1)+pos(2)*pos(2));
  float vertex_r_true = sqrt(posv(0)*posv(0)+posv(1)*posv(1)+posv(2)*posv(2));
  float diff_vertex_r = vertex_r_ring - vertex_r_true;
  m_hSvc.h1D("diff_vertex_r","","")->Fill(diff_vertex_r);
  m_hSvc.h1D(Form("diff_vertex_r_nring%d_mulike%d",nRing,n_mulike_angle),"","")->Fill(diff_vertex_r);

  if(kDebugMode){
    cout << "d_wall true/reco=" << wallv(0) << "/" << wall(0) << endl;
    cout << "true vertex position x/y/z=" << posv(0) << "/" << posv(1) << "/" << posv(2) << endl;
    cout << "reco vertex position x/y/z=" << pos(0) << "/" << pos(1) << "/" << pos(2) << endl;
    cout << "nring=" << nRing << endl;
    cout << "n_elike ip/pattern/angle=" << n_elike << "/" << n_elike_pattern << "/" << n_elike_angle << endl;
    cout << "n_mulike ip/pattern/angle=" << n_mulike << "/" << n_mulike_pattern << "/" << n_mulike_angle << endl;
    cout << "npar/npar2=" << npar(0) << "/" << npar2(0) << endl;
  }

  float min_mom=99999999, mid_mom=99999999, max_mom=99999999;
  int min_mom_id=-1,mid_mom_id=-1,max_mom_id=-1;
  TLorentzVector total_true_vec,total_gamma_true_vec;
  for(int v=1;v<4;v++){
    float expected_opening_angle=0;
    if(ipv(v)==2 || ipv(v)==3) expected_opening_angle = CalcOpeningAngle(0,pmomv(v));
    if(ipv(v)==5 || ipv(v)==6) expected_opening_angle = CalcOpeningAngle(1,pmomv(v));
    if(kDebugMode){
      cout << "true lepton pid/mom=" << ipv(v) << "/" << pmomv(v) << endl;
      cout << "true direction x/y/z=" << dirv(v,0) << "/" << dirv(v,1) << "/" << dirv(v,2) << endl;
      cout << "expected opening angle is " << expected_opening_angle << endl;
    }
    if(pmomv(v)<min_mom) {
      max_mom=mid_mom;
      mid_mom=min_mom;
      min_mom=pmomv(v);
      max_mom_id=mid_mom_id;
      mid_mom_id=min_mom_id;
      min_mom_id=v;
    }
    else if(pmomv(v)<mid_mom) {
      max_mom=mid_mom;
      mid_mom=pmomv(v);
      max_mom_id=mid_mom_id;
      mid_mom_id=v;
    }
    else if(pmomv(v)<max_mom) {
      max_mom=pmomv(v);
      max_mom_id=v;
    }
    TLorentzVector this_lepton = GetTLorentzVectorVector(v);
    total_true_vec = total_true_vec + this_lepton;
    float closest_angle=9999;
    int closest_ring_id=-1;
    for(int r=0;r<nRing;r++){
      TLorentzVector ring = GetTLorentzVectorRing(r,0);//just for angle calculation
      float prob_angle = sqrt(fabs(probms(r,1)))-sqrt(fabs(probms(r,2)));
      float prob_pattern = prmslg(r,1) - prmslg(r,2);
      int ip_angle = (prob_angle<0)? 2 : 3 ;
      int ip_pattern = (prob_pattern<0)? 2 : 3;
      float angle_lep_ring = this_lepton.Angle(ring.Vect())*180./3.14159;
      if(kDebugMode) {
        cout << "ring ip/pattern/angle=" << ip(r) << "/" << ip_pattern << "/" << ip_angle << endl;
        cout << "opening angle ang/e/m=" << ang(r) << "/" << ange(r) << "/" << angm(r) << endl;
        cout << "ring dir x/y/z=" << dir(r,0) << "/" << dir(r,1) << "/" << dir(r,2) << endl;
        cout << "momentum amom/e/m=" << amom(r) << "/" << amome(r) << "/" << amomm(r) << endl;
      }
      m_hSvc.h1D(Form("true_angle_lepton_and_ring_nring%d",nRing),"","")->Fill(angle_lep_ring);
      if(angle_lep_ring<20) {
        if(kDebugMode) cout << "match!!" << endl;
        if(angle_lep_ring<closest_angle){
          closest_angle = angle_lep_ring;
          closest_ring_id = r;
        }
      }
    }
    if(kDebugMode) cout << "closest_ring_id=" << closest_ring_id << endl;
    m_hSvc.h1D("true_mom_lepton","","")->Fill(pmomv(v));
    m_hSvc.h1D(Form("true_mom_lepton_nring%d",nRing),"","")->Fill(pmomv(v));
    if(ipv(v)==5 || ipv(v)==6) {
      m_hSvc.h1D("true_mom_muon","","")->Fill(pmomv(v));
      m_hSvc.h1D(Form("true_mom_muon_nring%d",nRing),"","")->Fill(pmomv(v));
      m_hSvc.h1D(Form("true_mom_muon_nring%d_mulike%d",nRing,n_mulike_angle),"","")->Fill(pmomv(v));
    }
    if(ipv(v)==2 || ipv(v)==3)m_hSvc.h1D("true_mom_electron","","")->Fill(pmomv(v));
    if(closest_ring_id!=-1){//true & ring matched
      float opening_angle = (ip(closest_ring_id)==2)? ange(closest_ring_id) : angm(closest_ring_id);
      float diff_opening_angle = opening_angle-expected_opening_angle;
      float diff_opening_angle_muon = angm(closest_ring_id) - expected_opening_angle;
      float diff_opening_angle_electron = ange(closest_ring_id) - expected_opening_angle;
      if(kDebugMode){
        if(wallv(0)>200 && fabs(diff_opening_angle)>10) cout << "diff_opening _angle is worse!!" << endl;
        cout << "diff_opening_angle=" << diff_opening_angle << endl;
      }
      float mom_pattern = ( (prmslg(closest_ring_id,1) - prmslg(closest_ring_id,2)) < 0 )? amome(closest_ring_id) : amomm(closest_ring_id);
      float mom_angle = ( sqrt(fabs(probms(closest_ring_id,1)))-sqrt(fabs(probms(closest_ring_id,2))) < 0 )? amome(closest_ring_id) : amomm(closest_ring_id);
      float residual_emom = (amome(closest_ring_id)-pmomv(v))/pmomv(v);
      float residual_mmom = (amomm(closest_ring_id)-pmomv(v))/pmomv(v);
      float residual_mom_angle = (mom_angle-pmomv(v))/pmomv(v);
      float residual_mom_pattern = (mom_pattern-pmomv(v))/pmomv(v);
      float prob_angle = sqrt(fabs(probms(closest_ring_id,1)))-sqrt(fabs(probms(closest_ring_id,2)));
      float prob_pattern = prmslg(closest_ring_id,1) - prmslg(closest_ring_id,2);
      m_hSvc.h1D("expected_opening_angle","","")->Fill(expected_opening_angle);
      m_hSvc.h2D("true_mom_expected_opening_angle","","")->Fill(pmomv(v),expected_opening_angle);
      m_hSvc.h2D("true_mom_opening_angle","","")->Fill(pmomv(v),opening_angle);
      m_hSvc.h1D("opening_angle","","")->Fill(opening_angle);
      m_hSvc.h1D("residual_emom","","")->Fill(residual_emom);
      m_hSvc.h1D("residual_mmom","","")->Fill(residual_mmom);
      m_hSvc.h1D("diff_opening_angle","","")->Fill(diff_opening_angle);
      m_hSvc.h1D(Form("diff_opening_angle_ip%d",ip(closest_ring_id)),"","")->Fill(diff_opening_angle);
      m_hSvc.h1D("true_mom_lepton_match_ring","","")->Fill(pmomv(v));
      m_hSvc.h1D(Form("true_mom_lepton_match_ring_nring%d",nRing),"","")->Fill(pmomv(v));
      if(ipv(v)==5 || ipv(v)==6) m_hSvc.h1D(Form("true_mom_muon_match_ring_nring%d",nRing),"","")->Fill(pmomv(v));
      if(prob_angle<0) {
        m_hSvc.h1D("true_mom_lepton_match_ring_angle_elike","","")->Fill(pmomv(v));
        m_hSvc.h1D(Form("true_mom_lepton_match_ring_angle_elike_nring%d",nRing),"","")->Fill(pmomv(v));
        m_hSvc.h1D("diff_opening_angle_angle_elike","","")->Fill(diff_opening_angle);
      }
      else {
        m_hSvc.h1D("true_mom_lepton_match_ring_angle_mulike","","")->Fill(pmomv(v));
        m_hSvc.h1D(Form("true_mom_lepton_match_ring_angle_mulike_nring%d",nRing),"","")->Fill(pmomv(v));
        if(ipv(v)==5 || ipv(v)==6) m_hSvc.h1D(Form("true_mom_muon_match_ring_angle_mulike_nring%d",nRing),"","")->Fill(pmomv(v));
        m_hSvc.h1D("diff_opening_angle_angle_mulike","","")->Fill(diff_opening_angle);
      }
      if(prob_pattern<0) {
        m_hSvc.h1D("true_mom_lepton_match_ring_charge_elike","","")->Fill(pmomv(v));
        m_hSvc.h1D("diff_opening_angle_charge_elike","","")->Fill(diff_opening_angle);
      }
      else {
        m_hSvc.h1D("true_mom_lepton_match_ring_charge_mulike","","")->Fill(pmomv(v));
        m_hSvc.h1D("diff_opening_angle_charge_mulike","","")->Fill(diff_opening_angle);
      }
      for(int p=0;p<6;p++){
        if(pmomv(v)>100*p && pmomv(v)<100+100*p){
          m_hSvc.h1D(Form("diff_opening_angle_mom%d_%d",100*p,100+100*p),"","")->Fill(diff_opening_angle);
          m_hSvc.h1D(Form("residual_emom_mom%d_%d",100*p,100+100*p),"","")->Fill(residual_emom);
          m_hSvc.h1D(Form("residual_mmom_mom%d_%d",100*p,100+100*p),"","")->Fill(residual_mmom);
          if(ipv(v)==2 || ipv(v)==3) {//electron
            m_hSvc.h1D(Form("diff_opening_angle_electron_mom%d_%d",100*p,100+100*p),"","")->Fill(diff_opening_angle_electron);
            m_hSvc.h1D(Form("prob_angle_electron_mom%d_%d",100*p,100+100*p),"","")->Fill(prob_angle);
          }
          if(ipv(v)==5 || ipv(v)==6) {//muon
            m_hSvc.h1D(Form("diff_opening_angle_muon_mom%d_%d",100*p,100+100*p),"","")->Fill(diff_opening_angle_muon);
            m_hSvc.h1D(Form("prob_angle_muon_mom%d_%d",100*p,100+100*p),"","")->Fill(prob_angle);
          }
        }
        if(diff_vertex_r>5*p && diff_vertex_r<5+5*p){
          m_hSvc.h1D(Form("diff_opening_angle_diff_vertex_r%d_%d",5*p,5+5*p),"","")->Fill(diff_opening_angle);
        }
      }
    }

    float closest_angle_lep_lep=9999;
    for(int vv=1;vv<4;vv++){
      if(vv==v) continue;
      TLorentzVector this_lepton2 = GetTLorentzVectorVector(vv);
      float angle_lep_lep = this_lepton.Angle(this_lepton2.Vect())*180./3.14159;
      //cout << "angle_lep_lep=" << angle_lep_lep << endl;
      if(closest_angle_lep_lep>angle_lep_lep) closest_angle_lep_lep=angle_lep_lep;
      if(pmomv(v)>150){
        m_hSvc.h1D("true_angle_lepton_and_lepton","","")->Fill(angle_lep_lep);
        if(closest_ring_id!=-1) m_hSvc.h1D("true_angle_lepton_and_lepton_match_ring","","")->Fill(angle_lep_lep);
      }
    }
    //cout << "closest_angle_lep_lep=" << closest_angle_lep_lep << endl;
    if(v==1){
      //cout << "Fill!!" << endl;
    }
  }
  TLorentzVector first_lep = GetTLorentzVectorVector(1);
  TLorentzVector second_lep = GetTLorentzVectorVector(2);
  TLorentzVector third_lep = GetTLorentzVectorVector(3);
  float mass_1st_2nd = (first_lep + second_lep).M();
  float mass_2nd_3rd = (second_lep + third_lep).M();
  float mass_3rd_1st = (third_lep + first_lep).M();
  float total_true_mass = total_true_vec.M();
  float total_true_mom = total_true_vec.P();
  float residual_total_mass = (total_mass-total_true_mass)/total_true_mass;
  float residual_total_mom = (total_mom-total_true_mom)/total_true_mom;
  m_hSvc.h1D(Form("total_true_mass_fp%d",is_free_proton),"","")->Fill(total_true_mass);
  m_hSvc.h1D(Form("total_true_mom_fp%d",is_free_proton),"","")->Fill(total_true_mom);
  m_hSvc.h2D(Form("total_true_mass_true_mom_fp%d",is_free_proton),"","")->Fill(total_true_mass,total_true_mom);
  m_hSvc.h2D(Form("total_true_mass_1st_lepton_mom_fp%d",is_free_proton),"","")->Fill(total_true_mass,pmomv(1));
  m_hSvc.h2D(Form("total_true_mom_1st_lepton_mom_fp%d",is_free_proton),"","")->Fill(total_true_mom,pmomv(1));
  m_hSvc.h1D("residual_total_mass","","")->Fill(residual_total_mass);
  m_hSvc.h1D("residual_total_mom","","")->Fill(residual_total_mom);
  m_hSvc.h1D(Form("residual_total_mass_nring%d_mulike%d",nRing,n_mulike_angle),"","")->Fill(residual_total_mass);
  m_hSvc.h1D(Form("residual_total_mom_nring%d_mulike%d",nRing,n_mulike_angle),"","")->Fill(residual_total_mom);
  m_hSvc.h1D(Form("true_mom_1st_lepton_fp%d",is_free_proton),"","")->Fill(pmomv(1));
  m_hSvc.h1D(Form("true_mom_2nd_lepton_fp%d",is_free_proton),"","")->Fill(pmomv(2));
  m_hSvc.h1D(Form("true_mom_3rd_lepton_fp%d",is_free_proton),"","")->Fill(pmomv(3));
  m_hSvc.h1D(Form("true_energy_1st_lepton_fp%d",is_free_proton),"","")->Fill(CalcEnergyVector(1));
  m_hSvc.h1D(Form("true_energy_2nd_lepton_fp%d",is_free_proton),"","")->Fill(CalcEnergyVector(2));
  m_hSvc.h1D(Form("true_energy_3rd_lepton_fp%d",is_free_proton),"","")->Fill(CalcEnergyVector(3));
  m_hSvc.h1D(Form("true_mass_1st_2nd_fp%d",is_free_proton),"","")->Fill(mass_1st_2nd);
  m_hSvc.h1D(Form("true_mass_2nd_3rd_fp%d",is_free_proton),"","")->Fill(mass_2nd_3rd);
  m_hSvc.h1D(Form("true_mass_3rd_1st_fp%d",is_free_proton),"","")->Fill(mass_3rd_1st);
  m_hSvc.graph(Form("true_mass_1st_2nd_vs_2nd_3rd_fp%d",is_free_proton),"","")->SetPoint(graph_point,mass_1st_2nd,mass_2nd_3rd);
  //cout << "total true vec mass/mom=" << total_true_vec.M() << "/" << total_true_vec.P() << endl;
  //cout << "final id min/mid/max=" << min_mom_id << "/" << mid_mom_id << "/" << max_mom_id << endl;
  //cout << "final mom min/mid/max=" << min_mom << "/" << mid_mom << "/" << max_mom << endl;
  TLorentzVector min_mom_lep = GetTLorentzVectorVector(min_mom_id);
  TLorentzVector mid_mom_lep = GetTLorentzVectorVector(mid_mom_id);
  TLorentzVector max_mom_lep = GetTLorentzVectorVector(max_mom_id);
  float angle_min_mid = min_mom_lep.Angle(mid_mom_lep.Vect())*180./3.14159;
  float angle_min_max = min_mom_lep.Angle(max_mom_lep.Vect())*180./3.14159;
  float angle_mid_max = mid_mom_lep.Angle(max_mom_lep.Vect())*180./3.14159;
  float min_angle=-1;
  if(angle_min_mid<angle_min_max){
    if(angle_mid_max<angle_min_mid) min_angle=angle_min_mid;
    else min_angle = angle_min_mid;
  }else{
    if(angle_mid_max<angle_min_max) min_angle=angle_mid_max;
    else min_angle = angle_min_max;
  }
  //cout << "min_angle=" << min_angle << endl;

  m_hSvc.h1D("true_min_mom_lepton","","")->Fill(min_mom);
  m_hSvc.h1D("true_mid_mom_lepton","","")->Fill(mid_mom);
  m_hSvc.h1D("true_max_mom_lepton","","")->Fill(max_mom);

  if(min_mom>200) m_hSvc.h1D(Form("true_min_angle_lepton_lepton_nring%d",nRing),"","")->Fill(min_angle);
  m_hSvc.h1D(Form("true_max_mom_lepton_nring%d",nRing),"","")->Fill(max_mom);
  m_hSvc.h1D(Form("true_mid_mom_lepton_nring%d",nRing),"","")->Fill(mid_mom);
  m_hSvc.h1D(Form("true_min_mom_lepton_nring%d",nRing),"","")->Fill(min_mom);
  m_hSvc.h1D(Form("true_angle_min_mid_lepton_nring%d",nRing),"","")->Fill(angle_min_mid);
  m_hSvc.h1D(Form("true_angle_min_max_lepton_nring%d",nRing),"","")->Fill(angle_min_max);
  m_hSvc.h1D(Form("true_angle_mid_max_lepton_nring%d",nRing),"","")->Fill(angle_mid_max);

  int n_true_decayE=0;
  total_gamma_true_vec = GetTLorentzVectorVector(1);
  for(int t=0;t<npar2(0);t++){
    if(kDebugMode) cout << "par2 pid/origin=" << ipv2(t) << "/" << iorg(t) << endl;
    if(ipv2(t)==2 || ipv2(t)==3) n_true_decayE++;
    if(ipv2(t)==1) {
      m_hSvc.h1D("true_mom_gamma","","")->Fill(pmomv2(t));
      total_gamma_true_vec = total_gamma_true_vec + GetTLorentzVectorVector2(t);
    }
  }
  m_hSvc.h1D(Form("nMulikeRing_angle_nring%d_trueDecayE%d",nRing,n_true_decayE),"","")->Fill(n_mulike_angle);
  float total_gamma_true_mass = total_gamma_true_vec.M();
  float total_gamma_true_mom = total_gamma_true_vec.P();
  float residual_total_gamma_mass = (total_mass-total_gamma_true_mass)/total_gamma_true_mass;
  float residual_total_gamma_mom = (total_mom-total_gamma_true_mom)/total_gamma_true_mom;
  m_hSvc.h1D("total_gamma_true_mass","","")->Fill(total_gamma_true_mass);
  m_hSvc.h1D("total_gamma_true_mom","","")->Fill(total_gamma_true_mom);
  m_hSvc.h1D(Form("total_gamma_true_mass_fp%d",is_free_proton),"","")->Fill(total_gamma_true_mass);
  m_hSvc.h1D(Form("total_gamma_true_mom_fp%d",is_free_proton),"","")->Fill(total_gamma_true_mom);
  m_hSvc.h1D("residual_total_gamma_mass","","")->Fill(residual_total_gamma_mass);
  m_hSvc.h1D("residual_total_gamma_mom","","")->Fill(residual_total_gamma_mom);
  m_hSvc.h1D(Form("residual_total_gamma_mass_nring%d_mulike%d",nRing,n_mulike_angle),"","")->Fill(residual_total_gamma_mass);
  m_hSvc.h1D(Form("residual_total_gamma_mom_nring%d_mulike%d",nRing,n_mulike_angle),"","")->Fill(residual_total_gamma_mom);

  if(kDebugMode) cout << "decayE true/detect=" << n_true_decayE << "/" << nDecayE << endl;
  m_hSvc.h1D("n_true_decayE","","")->Fill(n_true_decayE);
  m_hSvc.h1D(Form("nDecayE_true_decayE_%d",n_true_decayE),"","")->Fill(nDecayE);

  if(min_mom>200){
    expected_3ring_events_electron+=1.;
    expected_3ring_events_muon+=1.;
  }else{
    for(int p=0;p<10;p++){
      if(min_mom>p*20 && min_mom<20+p*20){
        expected_3ring_events_electron+=eff_e_ring[p];
      }
    }
  }
  //cout << "expected_3ring_events_electron=" << expected_3ring_events_electron << endl;
  graph_point++;//for tgraph
}

int OscNtupleManager::GetParType(int ring_id){
  float min_dr=99.;
  int min_pid=0;
  for(int p=1;p<npar(0);p++){//initial proton is skipped
    float dr = CalcDrRingVector(ring_id,p);
    //cout << "pid=" << ipv(p) << endl;
    //cout << "dr=" << dr << endl;
    if(/*dr<0.2 &&*/ dr<min_dr){
      min_dr = dr;
      min_pid = ipv(p);
    }
  }
  return min_pid;

}

bool OscNtupleManager::SetTrueMode( int &ans_mode){
  int abs_mode = abs(mode(0));
  int neutrino_type = -1;
  int this_mode = -1;
  if(ipnu(0)==12) neutrino_type=0;//nue
  if(ipnu(0)==-12) neutrino_type=1;//nuebar
  if(ipnu(0)==14) neutrino_type=2;//numu
  if(ipnu(0)==-14) neutrino_type=3;//numubar

  //CC Quasi-Elastic
  if(abs_mode==1 || abs_mode==2) this_mode=1;//nue n -> e- p 
  //CC single pi from delta resonance
  if(abs_mode==11) this_mode=2;//nue p -> e- p pi+ 
  if(abs_mode==12) this_mode=3;//nue n -> e- p pi0 
  if(abs_mode==13) this_mode=4;//nue n -> e- n pi+ 
  //CC coherent pi 
  if(abs_mode==16) this_mode=5;//nue O(16) -> e- pi+ O(16)
  //CC multiple pi 
  if(abs_mode==21) this_mode=6;//nue (p or n) -> e- (n or p) multi-pi
  //CC Eta from delta resonance
  if(abs_mode==22) this_mode=7;//nue n -> e- p Eta
  //CC K from delta resonance
  if(abs_mode==23) this_mode=8;//nue n -> e- lambda K+
  //CC DIS 
  if(abs_mode==26) this_mode=9;//nue (n or p) -> e- (n or p) meson
  //NC single pi from delta resonance
  if(abs_mode==31) this_mode=10;//nue n -> nue n pi0 
  if(abs_mode==32) this_mode=11;//nue p -> nue p pi0 
  if(abs_mode==33) this_mode=12;//nue n -> nue p pi- 
  if(abs_mode==34) this_mode=13;//nue p -> nue n pi+ 
  //NC coherent pi
  if(abs_mode==36) this_mode=14;//nue O(16) -> nue pi0 O(16)
  //NC elastic + gamma (?)
  if(abs_mode==38 || abs_mode==39) this_mode=15;//nue (n or p) -> nue (n or p) gamma
  //NC multiple pi
  if(abs_mode==41) this_mode=16;//nue (p or n) -> nue (n or p) multi-pi
  //NC Eta from delta resonance
  if(abs_mode==42) this_mode=17;//nue n -> nue n Eta
  if(abs_mode==43) this_mode=18;//nue p -> nue p Eta
  //NC K from delta resonance
  if(abs_mode==44) this_mode=19;//nue n -> nue lambda K0
  if(abs_mode==45) this_mode=20;//nue p -> nue lambda K+
  //NC DIS 
  if(abs_mode==46) this_mode=21;//nue (n or p) -> nue (n or p) meson
  //NC elastic
  if(abs_mode==51) this_mode=22;//nue p -> nue p 
  if(abs_mode==52) this_mode=23;//nue n -> nue n 

  ans_mode = 23*neutrino_type + this_mode;

  return true;
}


float OscNtupleManager::CalcOpeningAngle(int id, float mom){//0:electron 1:muon
  float n = 1.333;
  float mass=0;
  if(id==0) mass = 0.511;//MeV
  if(id==1) mass = 105.7;//MeV
  float theta = acos(sqrt(1+pow(mass,2)/pow(mom,2))/n);
  float angle = theta*180/3.14159;
  return angle;

}

float OscNtupleManager::CalcEnergyVector(int id){
  float mass=0.;
  if(ipv(id)==1) mass = 0;//gamma
  if(ipv(id)==2 || ipv(id)==3) mass = 0.511;//MeV
  if(ipv(id)==5 || ipv(id)==6) mass = 105.7;//MeV
  float energy = sqrt(pmomv(id)*pmomv(id) + mass*mass);
  return energy;
}

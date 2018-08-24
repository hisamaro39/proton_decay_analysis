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
      m_hSvc.create1D(Form("zenith_angle_type%d",i),"",10,-1,1);
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

  //basic plot
  m_hSvc.create1D("hist","hist",100,0,100);
  m_hSvc.create1D("ErmsHax","",100,0,5);
  m_hSvc.create1D("interaction_mode","",120,-60,60);

  m_hSvc.create1D("distance_to_wall_vector","",100,0,1000);

  m_hSvc.create1D("mass_proton_vector","",125,0,1250);
  m_hSvc.create1D("mom_proton_vector","",100,0,1000);
  m_hSvc.create1D("mom_first_proton_vector","",100,0,1000);
  m_hSvc.create2D("mass_mom_proton_vector","",125,0,1250,100,0,1000);

  m_hSvc.create1D("minMom_each_charged_leptons_vector","",40,0,800);
  m_hSvc.create1D("maxMom_each_charged_leptons_vector","",40,0,800);
  m_hSvc.create2D("minMom_maxMom_each_charged_leptons_vector","",40,0,800,40,0,800);

  m_hSvc.create1D("nElikeRing","",6,0,6);
  m_hSvc.create1D("nMulikeRing","",6,0,6);
  m_hSvc.create1D("dr_e_ring","",200,0,2);
  m_hSvc.create1D("dr_mu_ring","",200,0,2);
  m_hSvc.create1D("dr_2gamma_vector","",200,0,2);
  m_hSvc.create1D("dr_1st_e_2nd_e_vector","",200,0,2);
  m_hSvc.create1D("dr_1st_e_3rd_e_vector","",200,0,2);
  m_hSvc.create1D("dr_1st_e_gamma_vector","",20,0,2);
  m_hSvc.create1D("dr_2nd_e_gamma_vector","",20,0,2);
  m_hSvc.create1D("dr_3rd_e_gamma_vector","",20,0,2);
  m_hSvc.create1D("dr_first_mu_other_mu_vector","",200,0,2);
  m_hSvc.create1D("mom_pi0_vector","",100,0,1000);
  m_hSvc.create1D("mom_each_charged_leptons_vector","",25,0,500);
  m_hSvc.create1D("theta_each_charged_leptons_vector","",32,0,3.2);
  m_hSvc.create1D("phi_each_charged_leptons_vector","",64,-3.2,3.2);
  m_hSvc.create1D("mom_gamma_vector","",125,0,500);
  m_hSvc.create1D("mom_gamma_match_ring_vector","",25,0,500);
  m_hSvc.create1D("mom_gamma_vector_prob_ring_angle_elike_gamma_match_ring","",25,0,500);
  m_hSvc.create1D("mom_gamma_vector_prob_hit_elike_gamma_match_ring","",25,0,500);
  m_hSvc.create1D("mom_gamma_vector_pid_elike_gamma_match_ring","",25,0,500);
  m_hSvc.create1D("mom_gamma_vector_prob_ring_angle_mulike_gamma_match_ring","",25,0,500);
  m_hSvc.create1D("mom_gamma_vector_prob_hit_mulike_gamma_match_ring","",25,0,500);
  m_hSvc.create1D("mom_gamma_vector_pid_mulike_gamma_match_ring","",25,0,500);
  m_hSvc.create1D("mom_e_vector","",25,0,500);
  m_hSvc.create1D("mom_mu_vector","",25,0,500);
  m_hSvc.create1D("mom_e_match_ring_vector","",25,0,500);
  m_hSvc.create1D("mom_e_vector_prob_ring_angle_elike_e_match_ring","",25,0,500);
  m_hSvc.create1D("mom_e_vector_prob_hit_elike_e_match_ring","",25,0,500);
  m_hSvc.create1D("mom_e_vector_prob_ring_angle_mulike_e_match_ring","",25,0,500);
  m_hSvc.create1D("mom_e_vector_prob_hit_mulike_e_match_ring","",25,0,500);
  m_hSvc.create1D("mom_mu_match_ring_vector","",25,0,500);
  m_hSvc.create1D("mom_mu_vector_prob_ring_angle_elike_mu_match_ring","",25,0,500);
  m_hSvc.create1D("mom_mu_vector_prob_hit_elike_mu_match_ring","",25,0,500);
  m_hSvc.create1D("mom_mu_vector_prob_ring_angle_mulike_mu_match_ring","",25,0,500);
  m_hSvc.create1D("mom_mu_vector_prob_hit_mulike_mu_match_ring","",25,0,500);
  m_hSvc.create1D("mom_e_vector_pid_elike_e_match_ring","",25,0,500);
  m_hSvc.create1D("mom_e_vector_pid_mulike_e_match_ring","",25,0,500);
  m_hSvc.create1D("mom_mu_vector_pid_elike_mu_match_ring","",25,0,500);
  m_hSvc.create1D("mom_mu_vector_pid_mulike_mu_match_ring","",25,0,500);
  m_hSvc.create2D("mom_pi0_dr_2gamma_vector","",100,0,1000,200,0,2);
  int r_max=4,mu_max=4,p_max=4;
  if(process_mode=="p_epi" || process_mode=="p_mupi"){
   r_max=1;mu_max=1;p_max=1;
  }
  for(int c=0;c<10;c++){//cut type
    for(int r=0;r<r_max;r++){//cut pattern of nring
      for(int mu=0;mu<mu_max;mu++){//cut pattern of mu-like ring
        for(int p=0;p<p_max;p++){//cut pattern of michel electron
          m_hSvc.create1D(Form("distance_to_wall_cut%d_nring%d_mulike%d_michel%d",c,r,mu,p),"",100,0,1000);
          m_hSvc.create1D(Form("visible_energy_cut%d_nring%d_mulike%d_michel%d",c,r,mu,p),"",100,0,500);
          m_hSvc.create1D(Form("nhit_OD_cut%d_nring%d_mulike%d_michel%d",c,r,mu,p),"",100,0,100);
          m_hSvc.create1D(Form("nRing_cut%d_nring%d_mulike%d_michel%d",c,r,mu,p),"",10,0,10);
          m_hSvc.create1D(Form("n_michel_electron_cut%d_nring%d_mulike%d_michel%d",c,r,mu,p),"",5,0,5);
          m_hSvc.create1D(Form("distance_to_wall_weight_osc_cut%d_nring%d_mulike%d_michel%d",c,r,mu,p),"",100,0,1000);
          m_hSvc.create1D(Form("visible_energy_weight_osc_cut%d_nring%d_mulike%d_michel%d",c,r,mu,p),"",100,0,500);
          m_hSvc.create1D(Form("nhit_OD_weight_osc_cut%d_nring%d_mulike%d_michel%d",c,r,mu,p),"",100,0,100);
          m_hSvc.create1D(Form("nRing_weight_osc_cut%d_nring%d_mulike%d_michel%d",c,r,mu,p),"",10,0,10);
          m_hSvc.create1D(Form("n_michel_electron_weight_osc_cut%d_nring%d_mulike%d_michel%d",c,r,mu,p),"",5,0,5);
          m_hSvc.create1D(Form("mass_proton_reco_cut%d_nring%d_mulike%d_michel%d",c,r,mu,p),"",125,0,1250);
          m_hSvc.create1D(Form("mass_proton_reco_weight_mc_cut%d_nring%d_mulike%d_michel%d",c,r,mu,p),"",125,0,1250);
          m_hSvc.create1D(Form("mass_proton_reco_weight_osc_cut%d_nring%d_mulike%d_michel%d",c,r,mu,p),"",125,0,1250);
          m_hSvc.create1D(Form("mom_proton_reco_cut%d_nring%d_mulike%d_michel%d",c,r,mu,p),"",100,0,1000);
          m_hSvc.create1D(Form("mom_proton_reco_weight_mc_cut%d_nring%d_mulike%d_michel%d",c,r,mu,p),"",100,0,1000);
          m_hSvc.create1D(Form("mom_proton_reco_weight_osc_cut%d_nring%d_mulike%d_michel%d",c,r,mu,p),"",100,0,1000);
          m_hSvc.create1D(Form("mass_pi0_reco_cut%d_nring%d_mulike%d_michel%d",c,r,mu,p),"",50,0,500);
          m_hSvc.create1D(Form("mass_pi0_reco_weight_osc_cut%d_nring%d_mulike%d_michel%d",c,r,mu,p),"",50,0,500);
          m_hSvc.create1D(Form("mass_pi0_reco_log10_cut%d_nring%d_mulike%d_michel%d",c,r,mu,p),"",40,0,4);
          m_hSvc.create2D(Form("mass_mom_proton_reco_cut%d_nring%d_mulike%d_michel%d",c,r,mu,p),"",125,0,1250,100,0,1000);
          m_hSvc.create2D(Form("all_ring_mass_mom_proton_reco_cut%d_nring%d_mulike%d_michel%d",c,r,mu,p),"",125,0,1250,100,0,1000);
          m_hSvc.create2D(Form("all_mulike_mass_mom_proton_reco_cut%d_nring%d_mulike%d_michel%d",c,r,mu,p),"",125,0,1250,100,0,1000);
          m_hSvc.create2D(Form("mass_mom_proton_reco_weight_mc_cut%d_nring%d_mulike%d_michel%d",c,r,mu,p),"",125,0,1250,100,0,1000);
          m_hSvc.create2D(Form("mass_mom_proton_reco_weight_osc_cut%d_nring%d_mulike%d_michel%d",c,r,mu,p),"",125,0,1250,100,0,1000);
          for(int n=1;n<6;n++){
            m_hSvc.create1D(Form("nElikeRing_nring%d_cut%d_nring%d_mulike%d_michel%d",n,c,r,mu,p),"",6,0,6);
            m_hSvc.create1D(Form("nMulikeRing_nring%d_cut%d_nring%d_mulike%d_michel%d",n,c,r,mu,p),"",6,0,6);
            m_hSvc.create1D(Form("nElikeRing_nring%d_weight_osc_cut%d_nring%d_mulike%d_michel%d",n,c,r,mu,p),"",6,0,6);
            m_hSvc.create1D(Form("nMulikeRing_nring%d_weight_osc_cut%d_nring%d_mulike%d_michel%d",n,c,r,mu,p),"",6,0,6);
            m_hSvc.create1D(Form("mass_all_ring_reco_nring%d_cut%d_nring%d_mulike%d_michel%d",n,c,r,mu,p),"",125,0,1250);
            m_hSvc.create1D(Form("mom_all_ring_reco_nring%d_cut%d_nring%d_mulike%d_michel%d",n,c,r,mu,p),"",125,0,1250);
            m_hSvc.create1D(Form("mass_all_mulike_reco_nring%d_cut%d_nring%d_mulike%d_michel%d",n,c,r,mu,p),"",125,0,1250);
            m_hSvc.create1D(Form("mom_all_mulike_reco_nring%d_cut%d_nring%d_mulike%d_michel%d",n,c,r,mu,p),"",125,0,1250);
          }
          for(int e=1;e<4;e++){//# of e-like or mu-like ring
            m_hSvc.create1D(Form("mass_proton_reco_elike%d_cut%d_nring%d_mulike%d_michel%d",e,c,r,mu,p),"",125,0,1250);
            m_hSvc.create1D(Form("mass_proton_reco_mulike%d_cut%d_nring%d_mulike%d_michel%d",e,c,r,mu,p),"",125,0,1250);
            m_hSvc.create1D(Form("mom_proton_reco_elike%d_cut%d_nring%d_mulike%d_michel%d",e,c,r,mu,p),"",100,0,1000);
            m_hSvc.create1D(Form("mom_proton_reco_mulike%d_cut%d_nring%d_mulike%d_michel%d",e,c,r,mu,p),"",100,0,1000);
            m_hSvc.create2D(Form("mass_mom_proton_reco_elike%d_cut%d_nring%d_mulike%d_michel%d",e,c,r,mu,p),"",125,0,1250,100,0,1000);
            m_hSvc.create2D(Form("mass_mom_proton_reco_mulike%d_cut%d_nring%d_mulike%d_michel%d",e,c,r,mu,p),"",125,0,1250,100,0,1000);
            m_hSvc.create1D(Form("mass_pi0_reco_elike%d_cut%d_nring%d_mulike%d_michel%d",e,c,r,mu,p),"",50,0,500);
            m_hSvc.create1D(Form("mass_pi0_reco_elike%d_weight_osc_cut%d_nring%d_mulike%d_michel%d",e,c,r,mu,p),"",50,0,500);
            m_hSvc.create1D(Form("mass_pi0_reco_log10_elike%d_cut%d_nring%d_mulike%d_michel%d",e,c,r,mu,p),"",40,0,4);
            for(int r2=2;r2<4;r2++){
              m_hSvc.create1D(Form("mass_proton_reco_nring%d_elike%d_cut%d_nring%d_mulike%d_michel%d",r2,e,c,r,mu,p),"",125,0,1250);
              m_hSvc.create1D(Form("mass_proton_reco_nring%d_mulike%d_cut%d_nring%d_mulike%d_michel%d",r2,e,c,r,mu,p),"",125,0,1250);
            }
          }
        }
      }
    }
  }


  for(int n=1;n<6;n++){
    m_hSvc.create1D(Form("mom_pi0_vector_nring%d",n),"",100,0,1000);
    m_hSvc.create1D(Form("mom_gamma_vector_nring%d",n),"",100,0,1000);
    m_hSvc.create1D(Form("decay_typenring%d",n),"",7,0,7);
    m_hSvc.create1D(Form("minMom_each_charged_leptons_vector_nring%d",n),"",40,0,800);
    m_hSvc.create1D(Form("middleMom_each_charged_leptons_vector_nring%d",n),"",40,0,800);
    m_hSvc.create1D(Form("maxMom_each_charged_leptons_vector_nring%d",n),"",40,0,800);
    m_hSvc.create1D(Form("nElikeRing_nring%d",n),"",6,0,6);
    m_hSvc.create1D(Form("nMulikeRing_nring%d",n),"",6,0,6);
    m_hSvc.create1D(Form("diff_vertex_x_nring%d",n),"",100,-100,100);
    m_hSvc.create1D(Form("diff_vertex_y_nring%d",n),"",100,-100,100);
    m_hSvc.create1D(Form("diff_vertex_z_nring%d",n),"",100,-100,100);
    m_hSvc.create1D(Form("diff_vertex_r_nring%d",n),"",100,-100,100);
    for(int p=0;p<4;p++){
      m_hSvc.create1D(Form("ring_matched_pid_nring%d_elike%d",n,p),"",20,0,20);
    }
  }
  for(int ra=0;ra<n_range_match_gamma_mom-1;ra++){
    m_hSvc.create1D(Form("dtheta_gamma_mom%d_%d",range_match_gamma_mom[ra],range_match_gamma_mom[ra+1]),"",50,-0.2,0.2);
    m_hSvc.create1D(Form("dphi_gamma_mom%d_%d",range_match_gamma_mom[ra],range_match_gamma_mom[ra+1]),"",50,-0.4,0.4);
    m_hSvc.create1D(Form("dtheta_elike_gamma_mom%d_%d",range_match_gamma_mom[ra],range_match_gamma_mom[ra+1]),"",50,-0.2,0.2);
    m_hSvc.create1D(Form("dphi_elike_gamma_mom%d_%d",range_match_gamma_mom[ra],range_match_gamma_mom[ra+1]),"",50,-0.4,0.4);
    m_hSvc.create1D(Form("dtheta_mulike_gamma_mom%d_%d",range_match_gamma_mom[ra],range_match_gamma_mom[ra+1]),"",50,-0.2,0.2);
    m_hSvc.create1D(Form("dphi_mulike_gamma_mom%d_%d",range_match_gamma_mom[ra],range_match_gamma_mom[ra+1]),"",50,-0.4,0.4);
    m_hSvc.create1D(Form("dtheta_e_mom%d_%d",range_match_e_mom[ra],range_match_e_mom[ra+1]),"",50,-0.2,0.2);
    m_hSvc.create1D(Form("dphi_e_mom%d_%d",range_match_e_mom[ra],range_match_e_mom[ra+1]),"",50,-0.4,0.4);
    m_hSvc.create1D(Form("dtheta_elike_e_mom%d_%d",range_match_e_mom[ra],range_match_e_mom[ra+1]),"",50,-0.2,0.2);
    m_hSvc.create1D(Form("dphi_elike_e_mom%d_%d",range_match_e_mom[ra],range_match_e_mom[ra+1]),"",50,-0.4,0.4);
    m_hSvc.create1D(Form("dtheta_mulike_e_mom%d_%d",range_match_e_mom[ra],range_match_e_mom[ra+1]),"",50,-0.2,0.2);
    m_hSvc.create1D(Form("dphi_mulike_e_mom%d_%d",range_match_e_mom[ra],range_match_e_mom[ra+1]),"",50,-0.4,0.4);
    m_hSvc.create1D(Form("dtheta_mu_mom%d_%d",range_match_mu_mom[ra],range_match_mu_mom[ra+1]),"",50,-0.2,0.2);
    m_hSvc.create1D(Form("dphi_mu_mom%d_%d",range_match_mu_mom[ra],range_match_mu_mom[ra+1]),"",50,-0.4,0.4);
    m_hSvc.create1D(Form("dtheta_elike_mu_mom%d_%d",range_match_mu_mom[ra],range_match_mu_mom[ra+1]),"",50,-0.2,0.2);
    m_hSvc.create1D(Form("dphi_elike_mu_mom%d_%d",range_match_mu_mom[ra],range_match_mu_mom[ra+1]),"",50,-0.4,0.4);
    m_hSvc.create1D(Form("dtheta_mulike_mu_mom%d_%d",range_match_mu_mom[ra],range_match_mu_mom[ra+1]),"",50,-0.2,0.2);
    m_hSvc.create1D(Form("dphi_mulike_mu_mom%d_%d",range_match_mu_mom[ra],range_match_mu_mom[ra+1]),"",50,-0.4,0.4);
    m_hSvc.create1D(Form("residual_mom_gamma_mom%d_%d",range_match_gamma_mom[ra],range_match_gamma_mom[ra+1]),"",100,-1,1);
    m_hSvc.create1D(Form("residual_elike_mom_gamma_mom%d_%d",range_match_gamma_mom[ra],range_match_gamma_mom[ra+1]),"",100,-1,1);
    m_hSvc.create1D(Form("residual_mulike_mom_gamma_mom%d_%d",range_match_gamma_mom[ra],range_match_gamma_mom[ra+1]),"",100,-1,1);
    m_hSvc.create1D(Form("residual_mom_e_mom%d_%d",range_match_e_mom[ra],range_match_e_mom[ra+1]),"",100,-1,1);
    m_hSvc.create1D(Form("residual_elike_mom_e_mom%d_%d",range_match_e_mom[ra],range_match_e_mom[ra+1]),"",100,-1,1);
    m_hSvc.create1D(Form("residual_mulike_mom_e_mom%d_%d",range_match_e_mom[ra],range_match_e_mom[ra+1]),"",100,-1,1);
    m_hSvc.create1D(Form("residual_mom_mu_mom%d_%d",range_match_mu_mom[ra],range_match_mu_mom[ra+1]),"",100,-1,1);
    m_hSvc.create1D(Form("residual_elike_mom_mu_mom%d_%d",range_match_mu_mom[ra],range_match_mu_mom[ra+1]),"",100,-1,1);
    m_hSvc.create1D(Form("residual_mulike_mom_mu_mom%d_%d",range_match_mu_mom[ra],range_match_mu_mom[ra+1]),"",100,-1,1);
    m_hSvc.create1D(Form("prob_ring_angle_gamma_match_ringmom%d_%d",range_match_gamma_mom[ra],range_match_gamma_mom[ra+1]),"",60,-60,60);
    m_hSvc.create1D(Form("prob_hit_gamma_match_ringmom%d_%d",range_match_gamma_mom[ra],range_match_gamma_mom[ra+1]),"",50,-500,500);
    m_hSvc.create1D(Form("prob_ring_angle_e_match_ringmom%d_%d",range_match_e_mom[ra],range_match_e_mom[ra+1]),"",60,-60,60);
    m_hSvc.create1D(Form("prob_hit_e_match_ringmom%d_%d",range_match_e_mom[ra],range_match_e_mom[ra+1]),"",100,-500,500);
    m_hSvc.create1D(Form("prob_ring_angle_mu_match_ringmom%d_%d",range_match_mu_mom[ra],range_match_mu_mom[ra+1]),"",100,-20,20);
    m_hSvc.create1D(Form("prob_hit_mu_match_ringmom%d_%d",range_match_mu_mom[ra],range_match_mu_mom[ra+1]),"",100,-200,200);
  }

  for(int p=0;p<2;p++){
    m_hSvc.create1D(Form("energy_first_proton_%s_vector",proton_type[p].c_str()),"",125,0,1250);
    m_hSvc.create1D(Form("mass_e_pi0_%s_vector",proton_type[p].c_str()),"",125,0,1250);
    m_hSvc.create1D(Form("mom_e_pi0_%s_vector",proton_type[p].c_str()),"",100,0,1000);
    m_hSvc.create1D(Form("mass_e_pi0_zoom_%s_vector",proton_type[p].c_str()),"",200,800,1000);
    m_hSvc.create1D(Form("reco_mass_proton_%s_vector",proton_type[p].c_str()),"",125,0,1250);
    m_hSvc.create1D(Form("reco_mass_proton_wo_nucleus_%s_vector",proton_type[p].c_str()),"",125,0,1250);
    m_hSvc.create1D(Form("reco_mass_eee_%s_vector",proton_type[p].c_str()),"",125,0,1250);
    m_hSvc.create1D(Form("reco_mom_eee_%s_vector",proton_type[p].c_str()),"",100,0,1000);

    m_hSvc.create1D(Form("mom_pi0_vector%s",proton_type[p].c_str()),"",100,0,1000);
    m_hSvc.create1D(Form("mom_e_vector%s",proton_type[p].c_str()),"",100,0,1000);
    for(int t=0;t<7;t++){//try 7 types of p -> e+ e+ e-
      m_hSvc.create1D(Form("reco_mass_eee_type%d_%s_vector",t,proton_type[p].c_str()),"",125,0,1250);
    }
    m_hSvc.create1D(Form("mass_charged_leptons_vector_%s",proton_type[p].c_str()),"",125,0,1250);
    m_hSvc.create1D(Form("mom_charged_leptons_vector_%s",proton_type[p].c_str()),"",100,0,1000);
    m_hSvc.create1D(Form("maxMom_each_charged_leptons_vector_%s",proton_type[p].c_str()),"",40,0,800);
    m_hSvc.create1D(Form("minMom_each_charged_leptons_vector_%s",proton_type[p].c_str()),"",40,0,800);
  }

  m_hSvc.create1D("cut_flow","",10,0,10);
  m_hSvc.create1D("cut_flow_weight_osc","",10,0,10);
  for(int r=0;r<r_max;r++){
    for(int mu=0;mu<mu_max;mu++){
      for(int m=0;m<p_max;m++){
        m_hSvc.create1D(Form("cut_flow_nring%d_mulike%d_michel%d",r,mu,m),"",10,0,10);
      }
    }
  }


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
    interaction_type = -1;
    if(process_input=="fcmc" || process_input=="fcdt"){
      bool sc1 = SetEventType(process_input,interaction_type);
      if(process_input=="fcmc") {
        weight = live_time_weight;
        mc_weight = GetMCweight(); 
        osc_weight = Oscillate();
        weight *= weight*osc_weight;//osc_weight include mc_weight
      }
      if(process_mode=="fcmc"){
        MakeOscillationPlot();
        continue;
      }
    }

    if(kDebugMode) {
      std::cout << "weight mc is " << mc_weight << std::endl;
      std::cout << "weight osc is " << osc_weight << std::endl;
      std::cout << "weight=" << weight << std::endl;
      std::cout << "mode=" << mode(0) << std::endl;
      std::cout << "ipnu=" << ipnu(0) << std::endl;
      std::cout << "npar/numnu=" << npar(0) << "/" << numnu(0) << std::endl;
    }
    m_hSvc.h1D("interaction_mode","","")->Fill(mode(0));

    if(process_mode=="p_epi") Process_pepi();
    if(process_mode=="p_mupi") Process_pmupi();
    if(process_mode=="p_eee") Process_peee();
    if(process_mode=="p_mumumu") Process_pmumumu();
    if(process_mode=="p_emumu") Process_pemumu();
    if(process_mode=="p_muee") Process_pmuee();
  }

  std::cout << " Classification: " << std::endl
    << "   Written: " << Good << std::endl
    << "   Skipped: " << Bad << std::endl
    << std::endl;

  std::cout << "Event type of " << process_mode << std::endl;
  if(process_mode=="p_eee"){
    std::cout << "p->e+ e+ e-: " << n_eee << std::endl;
    std::cout << "p->e+ e+ e- gamma: " << n_eeeg << std::endl;
    std::cout << "p->e+ e+ e- proton: " << n_eeep << std::endl;
    std::cout << "p->e+ e+ e- neutron: " << n_eeen << std::endl;
    std::cout << "p->e+ e+ e- gamma proton: " << n_eeegp << std::endl;
    std::cout << "p->e+ e+ e- neutron proton: " << n_eeenp << std::endl;
    std::cout << "p->e+ e+ e- gamma neutron: " << n_eeegn << std::endl;
    int sum_bound = n_eee + n_eeeg + n_eeep + n_eeen + n_eeegp+n_eeenp+n_eeegn;
    std::cout << "bound proton events: " << sum_bound << std::endl; 
    std::cout << "free proton events: " << n_free_proton << std::endl;
    std::cout << "all sum is " << n_free_proton + sum_bound << std::endl;
  }
  if(process_mode=="p_epi"){
    std::cout << "no interaction: " << n_noint << std::endl;
    std::cout << "absorption: " << n_abs << std::endl;
    std::cout << "scattering: " << n_scat << std::endl;
    std::cout << "charge exchange: " << n_charge << std::endl;
    std::cout << "pi production: " << n_prod << std::endl;
    int sum_bound = n_noint + n_abs + n_scat + n_charge + n_prod;
    std::cout << "bound proton events: " << sum_bound << std::endl; 
    std::cout << "free proton events: " << n_free_proton << std::endl;
    std::cout << "all sum is " << n_free_proton + sum_bound << std::endl; 
    std::cout << "informatcion about truth & ring matching" << std::endl;
    std::cout << "e match events: " << n_match_e << std::endl;
    std::cout << "2 gamma match events: " << n_match_2gamma << std::endl;
    std::cout << "1 gamma match events: " << n_match_1gamma << std::endl;
    std::cout << "0 gamma match events: " << n_match_0gamma << std::endl;
  }
  if(process_mode=="p_mupi"){
    std::cout << "no interaction: " << n_noint << std::endl;
    std::cout << "absorption: " << n_abs << std::endl;
    std::cout << "scattering: " << n_scat << std::endl;
    std::cout << "charge exchange: " << n_charge << std::endl;
    std::cout << "pi production: " << n_prod << std::endl;
    int sum_bound = n_noint + n_abs + n_scat + n_charge + n_prod;
    std::cout << "bound proton events: " << sum_bound << std::endl; 
    std::cout << "free proton events: " << n_free_proton << std::endl;
    std::cout << "all sum is " << n_free_proton + sum_bound << std::endl; 
    std::cout << "informatcion about truth & ring matching" << std::endl;
    std::cout << "mu match events: " << n_match_mu << std::endl;
  }

  if(kMakeNtuple){
    TFile *file = new TFile(output_ntuple.c_str(),"recreate");
    otree->Write();
    file->Close();
  }
  else m_hSvc.WriteOutput();

}


void OscNtupleManager::Process_pepi(){

  TLorentzVector vec_e, vec_pi0, vec_g1, vec_g2, sum_decay_particles, sum_decay_particles_wo_nucleus;
  vector<TLorentzVector> decay_particles, decay_particles_wo_nucleus;
  decay_particles.clear();decay_particles_wo_nucleus.clear();
  int n_true_pi0=0,n_true_piplus_piminus=0,n_true_proton=0;
  for(int v=0;v<nPar;v++){
    int pidgen = ipv(v);
    if(pidgen==7) {
      vec_pi0 = GetTLorentzVectorVector(v);
      n_true_pi0++;
    }
    if(pidgen==2) vec_e = GetTLorentzVectorVector(v);
    if(pidgen==8 || pidgen==9) n_true_piplus_piminus++;
    if(pidgen==14) n_true_proton++;
    if(v==0){//first proton
      m_hSvc.h1D("mom_first_proton_vector","","")->Fill(pmomv(v));
    }else{//decay prticles
      decay_particles.push_back(GetTLorentzVectorVector(v));
      if(pidgen!=14 && pidgen!=13) decay_particles_wo_nucleus.push_back(GetTLorentzVectorVector(v));
    }
  }
  if(is_free_proton) {
    n_free_proton++;
    event_type=0;
  }
  else{
    if(nPar==3 && n_true_pi0==1) {//no interatcion
      n_noint++;
      event_type=1;
    }
    if(n_true_pi0==0 && n_true_piplus_piminus==0) {//absorption
      n_abs++;
      event_type=2;
    }
    if(n_true_pi0==0 && n_true_piplus_piminus==1) {//charge exchange
      n_charge++;
      event_type=3;
    }
    if(n_true_piplus_piminus + n_true_pi0 > 1) {//pi production
      n_prod++;
      event_type=4;
    }
    if(nPar>3 && n_true_pi0==1 && n_true_piplus_piminus==0) {//scattering
      n_scat++;
      event_type=5;
    }
  }
  float mass_epi0_vector = (vec_e+vec_pi0).M();
  float mom_epi0_vector = (vec_e+vec_pi0).P();
  if(event_type==0 || event_type==1 || event_type==5)
    m_hSvc.h1D(Form("mass_epi0_vector_%s",proton_type[is_free_proton].c_str()),"","")->Fill(mass_epi0_vector);
  m_hSvc.h1D("mom_epi0_vector","","")->Fill(mom_epi0_vector);
  m_hSvc.h2D("mass_mom_epi0_vector","","")->Fill(mass_epi0_vector,mom_epi0_vector);
  if(npar(0)==3 && ipv(1)==2 && ipv(2)==7){ 
    m_hSvc.h1D(Form("mass_e_pi0_%s_vector",proton_type[is_free_proton].c_str()),"","")->Fill(mass_epi0_vector);
    m_hSvc.h1D(Form("mass_e_pi0_zoom_%s_vector",proton_type[is_free_proton].c_str()),"","")->Fill(mass_epi0_vector);
    m_hSvc.h1D(Form("mom_e_pi0_%s_vector",proton_type[is_free_proton].c_str()),"","")->Fill(mom_epi0_vector);
  }
  for(unsigned int s=0;s<decay_particles.size();s++) sum_decay_particles += decay_particles[s];
  for(unsigned int s=0;s<decay_particles_wo_nucleus.size();s++) sum_decay_particles_wo_nucleus += decay_particles_wo_nucleus[s];
  m_hSvc.h1D(Form("reco_mass_proton_%s_vector",proton_type[is_free_proton].c_str()),"","")->Fill(sum_decay_particles.M());
  m_hSvc.h1D(Form("reco_mass_proton_wo_nucleus_%s_vector",proton_type[is_free_proton].c_str()),"","")->Fill(sum_decay_particles_wo_nucleus.M());
  m_hSvc.h1D(Form("energy_first_proton_%s_vector",proton_type[is_free_proton].c_str()),"","")->Fill(sqrt(938.3*938.3+pmomv(0)*pmomv(0)));

  if(event_type==0 || event_type==1){//only e+ and pi0  in final state
    if(kDebugMode){
      std::cout << "only e+ and pi0 in final state!!" << std::endl;
      std::cout << "nPar=" << nPar << std::endl;
      std::cout << "nRing=" << nRing << std::endl;
      std::cout << "vertex of ring x/y/z = " << pos(0) << "/" << pos(1) << "/" << pos(2) << std::endl;
      std::cout << "vertex of vector  x/y/z = " << posv(0) << "/" << posv(1) << "/" << posv(2) << std::endl;
    }

    vec_g1 = GetTLorentzVectorVector2(0);
    vec_g2 = GetTLorentzVectorVector2(1);
    TLorentzVector vec_g1_from_pi0 = vec_g1 - vec_pi0;
    TLorentzVector vec_g2_from_pi0 = vec_g2 - vec_pi0;
    if(kDebugMode){
      std::cout << "pi0 mom/theta/phi=" 
        << vec_pi0.P() << "/" << vec_pi0.Theta() << "/" << vec_pi0.Phi() << std::endl;
      std::cout << "g1 mom/theta/phi=" 
        << vec_g1.P() << "/" << vec_g1.Theta() << "/" << vec_g1.Phi() << std::endl;
      std::cout << "g2 mom/theta/phi=" 
        << vec_g2.P() << "/" << vec_g2.Theta() << "/" << vec_g2.Phi() << std::endl;
      std::cout << "g1_from_pi0 mom/theta/phi=" 
        << vec_g1_from_pi0.P() << "/" << vec_g1_from_pi0.Theta() << "/" << vec_g1_from_pi0.Phi() << std::endl;
      std::cout << "g2_from_pi0 mom/theta/phi=" 
        << vec_g2_from_pi0.P() << "/" << vec_g2_from_pi0.Theta() << "/" << vec_g2_from_pi0.Phi() << std::endl;
    }

    float dx_vertex = pos(0) - posv(0);
    float dy_vertex = pos(1) - posv(1);
    float dz_vertex = pos(2) - posv(2);
    float vertex_r_ring = sqrt(pos(0)*pos(0)+pos(1)*pos(1)+pos(2)*pos(2));
    float vertex_r_vector = sqrt(posv(0)*posv(0)+posv(1)*posv(1)+posv(2)*posv(2));
    float dr_vertex = vertex_r_ring - vertex_r_vector;
    m_hSvc.h1D(Form("diff_vertex_x_nring%d",nRing),"","")->Fill(dx_vertex);
    m_hSvc.h1D(Form("diff_vertex_y_nring%d",nRing),"","")->Fill(dx_vertex);
    m_hSvc.h1D(Form("diff_vertex_z_nring%d",nRing),"","")->Fill(dx_vertex);
    m_hSvc.h1D(Form("diff_vertex_r_nring%d",nRing),"","")->Fill(dx_vertex);
    for(int r=0;r<nRing;r++){
      float dr_e = CalcDrRingVector(r,1);
      float dr_g1 = CalcDrRingVector2(r,0);
      float dr_g2 = CalcDrRingVector2(r,1);
      m_hSvc.h1D("dr_e_ring","","")->Fill(dr_e);
      float prob_ring_angle = probms(r,2)-probms(r,3);
      float prob_hit = prmslg(r,1) - prmslg(r,2);
      if(kDebugMode){
        std::cout << "ring x/y/z=" << dir(r,0) << "/" <<  dir(r,1) << "/" << dir(r,2) << std::endl;
        std::cout << "ring e-like x/y/z=" << msdir(r,0,0) << "/" << msdir(r,1,0) << "/" << msdir(r,2,0) << std::endl;
        std::cout << "dr e/g1/g2=" << dr_e << "/" << dr_g1 << "/" << dr_g2 << std::endl;
        std::cout << "pid mine/ip=" << prob_ring_angle << "/" << ip(r) << std::endl;
      }
      if(dr_e<0.2) {
        if(kDebugMode){
          std::cout << "match positron!!" << std::endl;
          std::cout << "probms(2,r)/probms(3,r)=" << probms(2,r) << "/" << probms(3,r) << std::endl;
          std::cout << "prob_ring_angle=" << prob_ring_angle << std::endl;
          std::cout << "prmslg(2,r)/prmslg(3,r)=" << prmslg(2,r) << "/" << prmslg(3,r) << std::endl;
          std::cout << "prob_hit=" << prob_hit << std::endl;
        }
        m_hSvc.h1D("mom_e_match_ring_vector","","")->Fill(pmomv(1));
        if(prob_ring_angle<0) m_hSvc.h1D("mom_e_vector_prob_ring_angle_elike_e_match_ring","","")->Fill(pmomv(1));
        else m_hSvc.h1D("mom_e_vector_prob_ring_angle_mulike_e_match_ring","","")->Fill(pmomv(1));
        if(prob_hit<0) m_hSvc.h1D("mom_e_vector_prob_hit_elike_e_match_ring","","")->Fill(pmomv(1));
        else m_hSvc.h1D("mom_e_vector_prob_hit_mulike_e_match_ring","","")->Fill(pmomv(1));
        n_match_e++;
        TVector3 vector_e = GetTVectorVector(1);
        TVector3 vector_ring = GetTVectorRing(r);
        TVector3 vector_ring_elike = GetTVectorRingPID(r,1);
        TVector3 vector_ring_mulike = GetTVectorRingPID(r,2);
        float dtheta = vector_e.Theta() - vector_ring.Theta();
        float dphi = vector_e.DeltaPhi(vector_ring);
        float dtheta_elike = vector_e.Theta() - vector_ring_elike.Theta();
        float dphi_elike = vector_e.DeltaPhi(vector_ring_elike);
        float dtheta_mulike = vector_e.Theta() - vector_ring_mulike.Theta();
        float dphi_mulike = vector_e.DeltaPhi(vector_ring_mulike);
        float amom_residual = (pmomv(1)-amom(r))/pmomv(1);
        float amome_residual = (pmomv(1)-amome(r))/pmomv(1);
        float amomm_residual = (pmomv(1)-amomm(r))/pmomv(1);
        for(int ra=0;ra<n_range_match_e_mom-1;ra++){
          if(pmomv(1)>range_match_e_mom[ra] && pmomv(1)<range_match_e_mom[ra+1]){
            m_hSvc.h1D(Form("dtheta_e_mom%d_%d",range_match_e_mom[ra],range_match_e_mom[ra+1]),"","")->Fill(dtheta);
            m_hSvc.h1D(Form("dphi_e_mom%d_%d",range_match_e_mom[ra],range_match_e_mom[ra+1]),"","")->Fill(dphi);
            m_hSvc.h1D(Form("dtheta_elike_e_mom%d_%d",range_match_e_mom[ra],range_match_e_mom[ra+1]),"","")->Fill(dtheta_elike);
            m_hSvc.h1D(Form("dphi_elike_e_mom%d_%d",range_match_e_mom[ra],range_match_e_mom[ra+1]),"","")->Fill(dphi_elike);
            m_hSvc.h1D(Form("dtheta_mulike_e_mom%d_%d",range_match_e_mom[ra],range_match_e_mom[ra+1]),"","")->Fill(dtheta_mulike);
            m_hSvc.h1D(Form("dphi_mulike_e_mom%d_%d",range_match_e_mom[ra],range_match_e_mom[ra+1]),"","")->Fill(dphi_mulike);
            m_hSvc.h1D(Form("residual_mom_e_mom%d_%d",range_match_e_mom[ra],range_match_e_mom[ra+1]),"","")->Fill(amom_residual);
            m_hSvc.h1D(Form("residual_elike_mom_e_mom%d_%d",range_match_e_mom[ra],range_match_e_mom[ra+1]),"","")->Fill(amome_residual);
            m_hSvc.h1D(Form("residual_mulike_mom_e_mom%d_%d",range_match_e_mom[ra],range_match_e_mom[ra+1]),"","")->Fill(amomm_residual);
          }
        }
      }
    }
    m_hSvc.h1D("mom_e_vector","","")->Fill(pmomv(1));
    bool match_g1=false,match_g2=false;
    for(int p2=0;p2<npar2(0);p2++){
      m_hSvc.h1D("mom_gamma_vector","","")->Fill(pmomv2(p2));
      m_hSvc.h1D(Form("mom_gamma_vector_nring%d",nRing),"","")->Fill(pmomv2(p2));
      for(int r=0;r<nRing;r++){
        float dr_g = CalcDrRingVector2(r,p2);
        float prob_ring_angle = probms(r,2)-probms(r,3);
        float prob_hit = prmslg(r,1) - prmslg(r,2);
        if(dr_g<0.2){//match true gamma and ring
          TVector3 vector_gamma = GetTVectorVector2(p2);
          TVector3 vector_ring = GetTVectorRing(r);
          TVector3 vector_ring_elike = GetTVectorRingPID(r,1);
          TVector3 vector_ring_mulike = GetTVectorRingPID(r,2);
          float dtheta = vector_gamma.Theta() - vector_ring.Theta();
          float dphi = vector_gamma.DeltaPhi(vector_ring);
          float dtheta_elike = vector_gamma.Theta() - vector_ring_elike.Theta();
          float dphi_elike = vector_gamma.DeltaPhi(vector_ring_elike);
          float dtheta_mulike = vector_gamma.Theta() - vector_ring_mulike.Theta();
          float dphi_mulike = vector_gamma.DeltaPhi(vector_ring_mulike);
          float amom_residual = (pmomv2(p2)-amom(r))/pmomv2(p2);
          float amome_residual = (pmomv2(p2)-amome(r))/pmomv2(p2);
          float amomm_residual = (pmomv2(p2)-amomm(r))/pmomv2(p2);
          m_hSvc.h1D("dtheta_gamma","","")->Fill(dtheta);
          m_hSvc.h1D("mom_gamma_match_ring_vector","","")->Fill(pmomv2(p2));
          if(prob_ring_angle<0) m_hSvc.h1D("mom_gamma_vector_prob_ring_angle_elike_gamma_match_ring","","")->Fill(pmomv2(p2));
          else m_hSvc.h1D("mom_gamma_vector_prob_ring_angle_mulike_gamma_match_ring","","")->Fill(pmomv2(p2));
          if(prob_hit<0) m_hSvc.h1D("mom_gamma_vector_prob_hit_elike_gamma_match_ring","","")->Fill(pmomv2(p2));
          else m_hSvc.h1D("mom_gamma_vector_prob_hit_mulike_gamma_match_ring","","")->Fill(pmomv2(p2));
          if(ip(r)==2) m_hSvc.h1D("mom_gamma_vector_pid_elike_gamma_match_ring","","")->Fill(pmomv2(p2));
          else m_hSvc.h1D("mom_gamma_vector_pid_mulike_gamma_match_ring","","")->Fill(pmomv2(p2));
          for(int ra=0;ra<n_range_match_gamma_mom-1;ra++){
            if(pmomv2(p2)>range_match_gamma_mom[ra] && pmomv2(p2)<range_match_gamma_mom[ra+1]){
              m_hSvc.h1D(Form("dtheta_gamma_mom%d_%d",range_match_gamma_mom[ra],range_match_gamma_mom[ra+1]),"","")->Fill(dtheta);
              m_hSvc.h1D(Form("dphi_gamma_mom%d_%d",range_match_gamma_mom[ra],range_match_gamma_mom[ra+1]),"","")->Fill(dphi);
              m_hSvc.h1D(Form("dtheta_elike_gamma_mom%d_%d",range_match_gamma_mom[ra],range_match_gamma_mom[ra+1]),"","")->Fill(dtheta_elike);
              m_hSvc.h1D(Form("dphi_elike_gamma_mom%d_%d",range_match_gamma_mom[ra],range_match_gamma_mom[ra+1]),"","")->Fill(dphi_elike);
              m_hSvc.h1D(Form("dtheta_mulike_gamma_mom%d_%d",range_match_gamma_mom[ra],range_match_gamma_mom[ra+1]),"","")->Fill(dtheta_mulike);
              m_hSvc.h1D(Form("dphi_mulike_gamma_mom%d_%d",range_match_gamma_mom[ra],range_match_gamma_mom[ra+1]),"","")->Fill(dphi_mulike);
              m_hSvc.h1D(Form("residual_mom_gamma_mom%d_%d",range_match_gamma_mom[ra],range_match_gamma_mom[ra+1]),"","")->Fill(amom_residual);
              m_hSvc.h1D(Form("residual_elike_mom_gamma_mom%d_%d",range_match_gamma_mom[ra],range_match_gamma_mom[ra+1]),"","")->Fill(amome_residual);
              m_hSvc.h1D(Form("residual_mulike_mom_gamma_mom%d_%d",range_match_gamma_mom[ra],range_match_gamma_mom[ra+1]),"","")->Fill(amomm_residual);
              m_hSvc.h1D(Form("prob_ring_angle_gamma_match_ringmom%d_%d",range_match_gamma_mom[ra],range_match_gamma_mom[ra+1]),"","")->Fill(prob_ring_angle);
              m_hSvc.h1D(Form("prob_hit_gamma_match_ringmom%d_%d",range_match_gamma_mom[ra],range_match_gamma_mom[ra+1]),"","")->Fill(prob_hit);
            }
          }
          if(p2==0) match_g1=true;
          else match_g2=true;
        }
      }
    }
    if(match_g1 && match_g2) n_match_2gamma++;
    else if(match_g1 || match_g2) n_match_1gamma++;
    else n_match_0gamma++;

    float dr_g1g2 = CalcDrVector2Vector2(0,1);
    m_hSvc.h1D("dr_2gamma_vector","","")->Fill(dr_g1g2);
    m_hSvc.h1D("mom_pi0_vector","","")->Fill(pmomv(2));
    m_hSvc.h1D(Form("mom_pi0_vector_nring%d",nRing),"","")->Fill(pmomv(2));
    m_hSvc.h2D("mom_pi0_dr_2gamma_vector","","")->Fill(pmomv(2),dr_g1g2);
    m_hSvc.h1D(Form("mom_pi0_vector%s",proton_type[is_free_proton].c_str()),"","")->Fill(pmomv(2));
    m_hSvc.h1D(Form("mom_e_vector%s",proton_type[is_free_proton].c_str()),"","")->Fill(pmomv(1));
  }

  //apply selection here

  if(wallv(0)<200) return;

  pass_cut[0][0]=true;
  if ( evis(0) > 30.0 && wall(0) > 200.0 && nhitac(0) <  nhitac_cut[skgen] ){//FC & FV
    pass_cut[1][0]=true;
  }

  if(nring(0)==2 || nring(0)==3){//# of cherenkov ring
    pass_cut[2][0]=true;
  }

  if( (nring(0)==2 && n_elike==2) || (nring(0)==3 && n_elike==3) ){ // # of e-like ring
    pass_cut[3][0]=true;
  }

  if(ndcy(0)==0){ // michel (decay) electron cut
    pass_cut[4][0]=true;
  }

  vector<TLorentzVector> egamma_cand;
  egamma_cand.clear();
  for(int r=0;r<nRing;r++)
    if(ip(r)==2) egamma_cand.push_back(GetTLorentzVectorRing(r,0));

  //std::cout << "n_elike=" << n_elike << std::endl;
  int gamma1_id=-1,gamma2_id=-1;
  float min_diff=99999999;
  if(egamma_cand.size()>1){
    for(unsigned int r1=0;r1<egamma_cand.size()-1;r1++){
      for(unsigned int r2=r1+1;r2<egamma_cand.size();r2++){
        //std::cout << "r1/r2=" << r1 << "/" << r2 << std::endl;
        float mass_pi0_reco = (egamma_cand[r1] + egamma_cand[r2]).M();
        //std::cout << "mass_pi0_reco=" << mass_pi0_reco << std::endl;
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
  if(kDebugMode){
    std::cout << "id gamma1/gamma2=" << gamma1_id << "/" << gamma2_id << std::endl;
    std::cout << "closest_mass_pi0_reco=" << closest_mass_pi0_reco << std::endl;
  }
  if(n_elike>=2)m_hSvc.h1D("mass_pi0_reco","","")->Fill(closest_mass_pi0_reco,weight);

  if(nRing==3 && closest_mass_pi0_reco>85 && closest_mass_pi0_reco<185) pass_cut[5][0]=true;
  if(nRing==2) pass_cut[5][0]=true;

  TLorentzVector total_vec;
  if(n_elike==3) total_vec = egamma_cand[0] + egamma_cand[1] + egamma_cand[2];
  if(n_elike==2) total_vec = egamma_cand[0] + egamma_cand[1];
  total_mass = total_vec.M();
  total_mom = total_vec.P();
  if(total_mass>800 && total_mass<1050 && total_mom<100){//total mass & low momentum
    pass_cut[6][0]=true;
  }
  if(total_mass>800 && total_mass<1050 && total_mom>100 && total_mom<250){//total mass & high momentum
    pass_cut[7][0]=true;
  }

  TLorentzVector all_ring_vec;
  for(int r=0;r<nRing;r++){
    TLorentzVector temp_vec = 
      (ip(r)==2)? GetTLorentzVectorRing(r,0) : GetTLorentzVectorRing(r,2);
    all_ring_vec += temp_vec;
  }
  all_ring_mass = all_ring_vec.M();
  all_ring_mom = all_ring_vec.P();
  if(kDebugMode) std::cout << "all ring mass/mom=" << all_ring_mass << "/" << all_ring_mom << std::endl;

  MakeCutFlowValidate();

}

void OscNtupleManager::Process_pmupi(){
  std::cout << "ndcy/GetDecayE=" << ndcy(0) << "/" << llm->GetDecayE() << std::endl;

  TLorentzVector vec_mu, vec_pi0;
  int n_true_pi0=0,n_true_piplus_piminus=0,n_true_proton=0;
  for(int v=0;v<nPar;v++){
    int pidgen = ipv(v);
    if(pidgen==7) {
      vec_pi0 = GetTLorentzVectorVector(v);
      n_true_pi0++;
    }
    if(pidgen==5) vec_mu = GetTLorentzVectorVector(v);
    if(pidgen==8 || pidgen==9) n_true_piplus_piminus++;
    if(pidgen==14) n_true_proton++;
  }
  if(is_free_proton) {
    n_free_proton++;
    event_type=0;
  }
  else{
    if(nPar==3 && n_true_pi0==1) {//no interatcion
      n_noint++;
      event_type=1;
    }
    if(n_true_pi0==0 && n_true_piplus_piminus==0) {//absorption
      n_abs++;
      event_type=2;
    }
    if(n_true_pi0==0 && n_true_piplus_piminus==1) {//charge exchange
      n_charge++;
      event_type=3;
    }
    if(n_true_piplus_piminus + n_true_pi0 > 1) {//pi production
      n_prod++;
      event_type=4;
    }
    if(nPar>3 && n_true_pi0==1 && n_true_piplus_piminus==0) {//scattering
      n_scat++;
      event_type=5;
    }
  }

  if(event_type==0 || event_type==1){//only mu+ and pi0  in final state
    if(kDebugMode){
      std::cout << "only e+ and pi0 in final state!!" << std::endl;
      std::cout << "nPar=" << nPar << std::endl;
      std::cout << "nRing=" << nRing << std::endl;
    }
    m_hSvc.h1D("mom_mu_vector","","")->Fill(pmomv(1));
    for(int r=0;r<nRing;r++){
      float dr_mu = CalcDrRingVector(r,1);
      m_hSvc.h1D("dr_mu_ring","","")->Fill(dr_mu);
      float prob_ring_angle = probms(r,2)-probms(r,3);
      float prob_hit = prmslg(r,1) - prmslg(r,2);
      if(dr_mu<0.2) {
        m_hSvc.h1D("prob_ring_angle_mu_match_ring","","")->Fill(prob_ring_angle);
        m_hSvc.h1D("prob_hit_mu_match_ring","","")->Fill(prob_hit);
        m_hSvc.h1D("mom_mu_match_ring_vector","","")->Fill(pmomv(1));
        if(prob_ring_angle<0) m_hSvc.h1D("mom_mu_vector_prob_ring_angle_elike_mu_match_ring","","")->Fill(pmomv(1));
        else m_hSvc.h1D("mom_mu_vector_prob_ring_angle_mulike_mu_match_ring","","")->Fill(pmomv(1));
        if(prob_hit<0) m_hSvc.h1D("mom_mu_vector_prob_hit_elike_mu_match_ring","","")->Fill(pmomv(1));
        else m_hSvc.h1D("mom_mu_vector_prob_hit_mulike_mu_match_ring","","")->Fill(pmomv(1));
        n_match_mu++;
        TVector3 vector_mu = GetTVectorVector(1);
        TVector3 vector_ring = GetTVectorRing(r);
        TVector3 vector_ring_elike = GetTVectorRingPID(r,1);
        TVector3 vector_ring_mulike = GetTVectorRingPID(r,2);
        float dtheta = vector_mu.Theta() - vector_ring.Theta();
        float dphi = vector_mu.DeltaPhi(vector_ring);
        float dtheta_elike = vector_mu.Theta() - vector_ring_elike.Theta();
        float dphi_elike = vector_mu.DeltaPhi(vector_ring_elike);
        float dtheta_mulike = vector_mu.Theta() - vector_ring_mulike.Theta();
        float dphi_mulike = vector_mu.DeltaPhi(vector_ring_mulike);
        float amom_residual = (pmomv(1)-amom(r))/pmomv(1);
        float amome_residual = (pmomv(1)-amome(r))/pmomv(1);
        float amomm_residual = (pmomv(1)-amomm(r))/pmomv(1);
        for(int ra=0;ra<n_range_match_mu_mom-1;ra++){
          if(pmomv(1)>range_match_mu_mom[ra] && pmomv(1)<range_match_mu_mom[ra+1]){
            m_hSvc.h1D(Form("dtheta_mu_mom%d_%d",range_match_mu_mom[ra],range_match_mu_mom[ra+1]),"","")->Fill(dtheta);
            m_hSvc.h1D(Form("dphi_mu_mom%d_%d",range_match_mu_mom[ra],range_match_mu_mom[ra+1]),"","")->Fill(dphi);
            m_hSvc.h1D(Form("dtheta_elike_mu_mom%d_%d",range_match_mu_mom[ra],range_match_mu_mom[ra+1]),"","")->Fill(dtheta_elike);
            m_hSvc.h1D(Form("dphi_elike_mu_mom%d_%d",range_match_mu_mom[ra],range_match_mu_mom[ra+1]),"","")->Fill(dphi_elike);
            m_hSvc.h1D(Form("dtheta_mulike_mu_mom%d_%d",range_match_mu_mom[ra],range_match_mu_mom[ra+1]),"","")->Fill(dtheta_mulike);
            m_hSvc.h1D(Form("dphi_mulike_mu_mom%d_%d",range_match_mu_mom[ra],range_match_mu_mom[ra+1]),"","")->Fill(dphi_mulike);
            m_hSvc.h1D(Form("residual_mom_mu_mom%d_%d",range_match_mu_mom[ra],range_match_mu_mom[ra+1]),"","")->Fill(amom_residual);
            m_hSvc.h1D(Form("residual_elike_mom_mu_mom%d_%d",range_match_mu_mom[ra],range_match_mu_mom[ra+1]),"","")->Fill(amome_residual);
            m_hSvc.h1D(Form("residual_mulike_mom_mu_mom%d_%d",range_match_mu_mom[ra],range_match_mu_mom[ra+1]),"","")->Fill(amomm_residual);
          }
        }
      }
    }
  }

  float mass_proton_vector = (vec_mu+vec_pi0).M();
  float mom_proton_vector = (vec_mu+vec_pi0).P();
  m_hSvc.h1D("mass_proton_vector","","")->Fill(mass_proton_vector);
  m_hSvc.h1D("mom_proton_vector","","")->Fill(mom_proton_vector);
  m_hSvc.h2D("mass_mom_proton_vector","","")->Fill(mass_proton_vector,mom_proton_vector);

  //apply selection here

  if(wallv(0)<200) return;

  pass_cut[0][0]=true;
  if ( evis(0) > 30.0 && wall(0) > 200.0 && nhitac(0) <  nhitac_cut[skgen] ){//FC & FV
    pass_cut[1][0]=true;
  }

  if(nring(0)==2 || nring(0)==3){//# of cherenkov ring
    pass_cut[2][0]=true;
  }

  if( n_mulike==1 ){ // # of mu-like ring
    pass_cut[3][0]=true;
  }

  if(ndcy(0)==1){ // michel (decay) electron cut
    pass_cut[4][0]=true;
  }

  vector<TLorentzVector> gamma_cand;
  TLorentzVector mu_cand;
  gamma_cand.clear();
  for(int r=0;r<nRing;r++){
    if(ip(r)==2) gamma_cand.push_back(GetTLorentzVectorRing(r,0));
    if(ip(r)==3) mu_cand = GetTLorentzVectorRing(r,2);
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
  if(kDebugMode){
    std::cout << "id gamma1/gamma2=" << gamma1_id << "/" << gamma2_id << std::endl;
    std::cout << "closest_mass_pi0_reco=" << closest_mass_pi0_reco << std::endl;
  }
  if(n_elike>=2) m_hSvc.h1D("mass_pi0_reco","","")->Fill(closest_mass_pi0_reco,weight);

  if(nRing==3 && closest_mass_pi0_reco>85 && closest_mass_pi0_reco<185) pass_cut[5][0]=true;
  if(nRing==2) pass_cut[5][0]=true;

  TLorentzVector total_vec;
  if(n_elike==2) total_vec = gamma_cand[0] + gamma_cand[1] + mu_cand;
  if(n_elike==1) total_vec = gamma_cand[0] + mu_cand;
  total_mass = total_vec.M();
  total_mom = total_vec.P();
  if(total_mass>800 && total_mass<1050 && total_mom<100){//total mass & low momentum
    pass_cut[6][0]=true;
  }
  if(total_mass>800 && total_mass<1050 && total_mom>100 && total_mom<250){//total mass & high momentum
    pass_cut[7][0]=true;
  }

  TLorentzVector all_ring_vec;
  for(int r=0;r<nRing;r++){
    TLorentzVector temp_vec = 
      (ip(r)==2)? GetTLorentzVectorRing(r,0) : GetTLorentzVectorRing(r,2);
    all_ring_vec += temp_vec;
  }
  all_ring_mass = all_ring_vec.M();
  all_ring_mom = all_ring_vec.P();
  if(kDebugMode) std::cout << "all ring mass/mom=" << all_ring_mass << "/" << all_ring_mom << std::endl;

  MakeCutFlowValidate();

}

void OscNtupleManager::Process_peee(){

  int n_true_e=0,n_true_g=0,n_true_p=0,n_true_n=0;
  float min_mom=99999999, min_mom_2nd=99999999, min_mom_3rd=99999999;
  for(int v=0;v<nPar;v++){
    int pidgen = ipv(v);
    if(pidgen==2 || pidgen==3) n_true_e++;
    if(pidgen==14) n_true_p++;
    if(pidgen==1) n_true_g++;
    if(pidgen==13) n_true_n++;
    if(v==1 || v==2 || v==3) {
      m_hSvc.h1D("mom_e_vector","","")->Fill(pmomv(v));
      m_hSvc.h1D("mom_each_charged_leptons_vector","","")->Fill(pmomv(v));
      TLorentzVector temp_lep = GetTLorentzVectorVector(v);
      m_hSvc.h1D("theta_each_charged_leptons_vector","","")->Fill(temp_lep.Theta());
      m_hSvc.h1D("phi_each_charged_leptons_vector","","")->Fill(temp_lep.Phi());
      if(pmomv(v)<min_mom) {
        min_mom_2nd=min_mom;
        min_mom=pmomv(v);
      }
      else if(pmomv(v)<min_mom_2nd) {
        min_mom_3rd=min_mom_2nd;
        min_mom_2nd=pmomv(v);
      }
      else if(pmomv(v)<min_mom_3rd) {
        min_mom_3rd=pmomv(v);
      }
    }
    for (int r=0;r<nRing;r++){
      float dr_e = CalcDrRingVector(r,v);
      float prob_ring_angle = probms(r,2)-probms(r,3);
      float prob_hit = prmslg(r,1) - prmslg(r,2);
      if(v==1 || v==2 || v==3){
        if(dr_e<0.2) {
          m_hSvc.h1D("mom_e_match_ring_vector","","")->Fill(pmomv(v));
          if(prob_ring_angle<0) m_hSvc.h1D("mom_e_vector_prob_ring_angle_elike_e_match_ring","","")->Fill(pmomv(v));
          else m_hSvc.h1D("mom_e_vector_prob_ring_angle_mulike_e_match_ring","","")->Fill(pmomv(v));
          if(prob_hit<0) m_hSvc.h1D("mom_e_vector_prob_hit_elike_e_match_ring","","")->Fill(pmomv(v));
          else m_hSvc.h1D("mom_e_vector_prob_hit_mulike_e_match_ring","","")->Fill(pmomv(v));
          if(ip(r)==2) m_hSvc.h1D("mom_e_vector_pid_elike_e_match_ring","","")->Fill(pmomv(v));
          else m_hSvc.h1D("mom_e_vector_pid_mulike_e_match_ring","","")->Fill(pmomv(v));
          for(int ra=0;ra<n_range_match_e_mom-1;ra++){
            if(pmomv(v)>range_match_e_mom[ra] && pmomv(v)<range_match_e_mom[ra+1]){
              m_hSvc.h1D(Form("prob_ring_angle_e_match_ringmom%d_%d",range_match_e_mom[ra],range_match_e_mom[ra+1]),"","")->Fill(prob_ring_angle);
              m_hSvc.h1D(Form("prob_hit_e_match_ringmom%d_%d",range_match_e_mom[ra],range_match_e_mom[ra+1]),"","")->Fill(prob_hit);
            }
          }
          TVector3 vector_e = GetTVectorVector(v);
          TVector3 vector_ring = GetTVectorRing(r);
          float dtheta = vector_e.Theta() - vector_ring.Theta();
          float dphi = vector_e.DeltaPhi(vector_ring);
          float amome_residual = (pmomv(v)-amome(r))/pmomv(v);
          for(int ra=0;ra<n_range_match_e_mom-1;ra++){
            if(pmomv(v)>range_match_e_mom[ra] && pmomv(v)<range_match_e_mom[ra+1]){
              m_hSvc.h1D(Form("dtheta_e_mom%d_%d",range_match_e_mom[ra],range_match_e_mom[ra+1]),"","")->Fill(dtheta);
              m_hSvc.h1D(Form("dphi_e_mom%d_%d",range_match_e_mom[ra],range_match_e_mom[ra+1]),"","")->Fill(dphi);
              m_hSvc.h1D(Form("residual_elike_mom_e_mom%d_%d",range_match_e_mom[ra],range_match_e_mom[ra+1]),"","")->Fill(amome_residual);
            }
          }
        }
      }
    }
  }
  if(kDebugMode) std::cout << "min mom 1st/2nd/3rd=" << min_mom << "/" << min_mom_2nd << "/" << min_mom_3rd << std::endl;
  m_hSvc.h1D("nElikeRing","","")->Fill(n_elike);
  m_hSvc.h1D(Form("nElikeRing_nring%d",nRing),"","")->Fill(n_elike);
  m_hSvc.h2D("minMom_maxMom_each_charged_leptons_vector","","")->Fill(min_mom, min_mom_3rd);
  m_hSvc.h1D("minMom_each_charged_leptons_vector","","")->Fill(min_mom);
  m_hSvc.h1D("maxMom_each_charged_leptons_vector","","")->Fill(min_mom_3rd);
  m_hSvc.h1D(Form("minMom_each_charged_leptons_vector_%s",proton_type[is_free_proton].c_str()),"","")->Fill(min_mom);
  m_hSvc.h1D(Form("maxMom_each_charged_leptons_vector_%s",proton_type[is_free_proton].c_str()),"","")->Fill(min_mom_3rd);
  m_hSvc.h1D(Form("minMom_each_charged_leptons_vector_nring%d",nRing),"","")->Fill(min_mom);
  m_hSvc.h1D(Form("middleMom_each_charged_leptons_vector_nring%d",nRing),"","")->Fill(min_mom_2nd);
  m_hSvc.h1D(Form("maxMom_each_charged_leptons_vector_nring%d",nRing),"","")->Fill(min_mom_3rd);

  if(is_free_proton) n_free_proton++;
  else{
    if(nPar==4) {
      n_eee++;
      event_type=0;
    }
    if(nPar==5 && n_true_g==1) {
      n_eeeg++;
      event_type=1;
    }
    if(nPar==5 && n_true_p==2) {
      n_eeep++;
      event_type=2;
    }
    if(nPar==5 && n_true_n==1) {
      n_eeen++;
      event_type=3;
    }
    if(nPar==6 && n_true_p==2 && n_true_g==1) {
      n_eeegp++;
      event_type=4;
    }
    if(nPar==6 && n_true_p==2 && n_true_n==1) {
      n_eeenp++;
      event_type=5;
    }
    if(nPar==6 && n_true_g==1 && n_true_n==1) {
      n_eeegn++;
      event_type=6;
    }
  }

  TLorentzVector vec_e1,vec_e2,vec_e3;
  vec_e1 = GetTLorentzVectorVector(1);
  vec_e2 = GetTLorentzVectorVector(2);
  vec_e3 = GetTLorentzVectorVector(3);
  float mass_eee_vector = (vec_e1+vec_e2+vec_e3).M();
  float mom_eee_vector = (vec_e1+vec_e2+vec_e3).P();
  m_hSvc.h1D(Form("mass_charged_leptons_vector_%s",proton_type[is_free_proton].c_str()),"","")->Fill(mass_eee_vector);
  m_hSvc.h1D(Form("mom_charged_leptons_vector_%s",proton_type[is_free_proton].c_str()),"","")->Fill(mom_eee_vector);
  m_hSvc.h1D("mom_eee_vector","","")->Fill(mom_eee_vector);
  m_hSvc.h2D("mass_mom_eee_vector","","")->Fill(mass_eee_vector,mom_eee_vector);
  m_hSvc.h1D(Form("reco_mass_eee_%s_vector",proton_type[is_free_proton].c_str()),"","")->Fill(mass_eee_vector);
  m_hSvc.h1D(Form("reco_mom_eee_%s_vector",proton_type[is_free_proton].c_str()),"","")->Fill(mom_eee_vector);
  m_hSvc.h1D(Form("reco_mass_eee_type%d_%s_vector",event_type,proton_type[is_free_proton].c_str()),"","")->Fill(mass_eee_vector);

  float dr_1st_e_2nd_e = CalcDrVectorVector(1,2);
  float dr_1st_e_3rd_e = CalcDrVectorVector(1,3);
  m_hSvc.h1D("nRing","","")->Fill(nring(0));
  m_hSvc.h1D(Form("decay_typenring%d",nRing),"","")->Fill(event_type);
  m_hSvc.h1D("dr_1st_e_2nd_e_vector","","")->Fill(dr_1st_e_2nd_e);
  m_hSvc.h1D("dr_1st_e_3rd_e_vector","","")->Fill(dr_1st_e_3rd_e);
  if(is_free_proton) m_hSvc.h1D(Form("decay_typenring%d",nRing),"","")->Fill(0);
  for(int nn=0;nn<nPar;nn++){
    if(ipv(nn)==1){
      float dr_1st_e_gamma = CalcDrVectorVector(1,nn);
      float dr_2nd_e_gamma = CalcDrVectorVector(2,nn);
      float dr_3rd_e_gamma = CalcDrVectorVector(3,nn);
      m_hSvc.h1D("dr_1st_e_gamma_vector","","")->Fill(dr_1st_e_gamma);
      m_hSvc.h1D("dr_2nd_e_gamma_vector","","")->Fill(dr_2nd_e_gamma);
      m_hSvc.h1D("dr_3rd_e_gamma_vector","","")->Fill(dr_3rd_e_gamma);
    }
  }

  for (int r=0;r<nRing;r++){
    int par_id = GetParType(r);
    std::cout << "matched pid = " << par_id << std::endl;
    m_hSvc.h1D(Form("ring_matched_pid_nring%d_elike%d",nRing,n_elike),"","")->Fill(par_id);
  }

  //apply selection here

  if(wallv(0)<200) return;

  pass_cut[0][0]=true;
  if ( evis(0) > 30.0 && wall(0) > 200.0 && nhitac(0) <  nhitac_cut[skgen] ){//FC & FV
    pass_cut[1][0]=true;
  }

  //# of cherenkov ring
  if(nring(0)==2) pass_cut[2][0]=true;
  if(nring(0)==3) pass_cut[2][1]=true;
  if(nring(0)==4) pass_cut[2][2]=true;

  //# of e-like ring
  if(n_mulike==0) pass_cut[3][0]=true;
  if(n_mulike==1) pass_cut[3][1]=true;
  if(n_mulike==2) pass_cut[3][2]=true;
  if(n_mulike==3) pass_cut[3][3]=true;

  // michel (decay) electron cut
  if(ndcy(0)==0) pass_cut[4][0]=true;
  if(ndcy(0)==1) pass_cut[4][1]=true;
  if(ndcy(0)==2) pass_cut[4][2]=true;
  if(ndcy(0)==3) pass_cut[4][3]=true;

  vector<TLorentzVector> e_cand;
  e_cand.clear();
  for(int r=0;r<nRing;r++)
    if(ip(r)==2) e_cand.push_back(GetTLorentzVectorRing(r,1));

  pass_cut[5][0] = true;//no pi0 selection
  TLorentzVector total_vec;
  if(n_elike==3) total_vec = e_cand[0] + e_cand[1] + e_cand[2];
  if(n_elike==2) total_vec = e_cand[0] + e_cand[1];
  total_mass = total_vec.M();
  total_mom = total_vec.P();
  if(total_mass>800 && total_mass<1050 && total_mom<100){//total mass & low momentum
    pass_cut[6][0]=true;
  }
  if(total_mass>800 && total_mass<1050 && total_mom>100 && total_mom<250){//total mass & high momentum
    pass_cut[7][0]=true;
  }

  TLorentzVector all_ring_vec;
  for(int r=0;r<nRing;r++){
    TLorentzVector temp_vec = 
      (ip(r)==2)? GetTLorentzVectorRing(r,0) : GetTLorentzVectorRing(r,2);
    all_ring_vec += temp_vec;
  }
  all_ring_mass = all_ring_vec.M();
  all_ring_mom = all_ring_vec.P();
  if(kDebugMode) std::cout << "all ring mass/mom=" << all_ring_mass << "/" << all_ring_mom << std::endl;

  MakeCutFlow();

}

void OscNtupleManager::Process_pmumumu(){

  float min_mom=99999999, min_mom_2nd=99999999, min_mom_3rd=99999999;
  for(int v=0;v<nPar;v++){
    int pidgen = ipv(v);
    if(v==1 || v==2 || v==3) {
      m_hSvc.h1D("mom_each_charged_leptons_vector","","")->Fill(pmomv(v));
      TLorentzVector temp_lep = GetTLorentzVectorVector(v);
      m_hSvc.h1D("theta_each_charged_leptons_vector","","")->Fill(temp_lep.Theta());
      m_hSvc.h1D("phi_each_charged_leptons_vector","","")->Fill(temp_lep.Phi());
      m_hSvc.h1D("mom_mu_vector","","")->Fill(pmomv(v));
      if(pmomv(v)<min_mom) {
        min_mom_2nd=min_mom;
        min_mom=pmomv(v);
      }
      else if(pmomv(v)<min_mom_2nd) {
        min_mom_3rd=min_mom_2nd;
        min_mom_2nd=pmomv(v);
      }
      else if(pmomv(v)<min_mom_3rd) {
        min_mom_3rd=pmomv(v);
      }
    }
    for (int r=0;r<nRing;r++){
      float dr_mu = CalcDrRingVector(r,v);
      float prob_ring_angle = probms(r,2)-probms(r,3);
      float prob_hit = prmslg(r,1) - prmslg(r,2);
      if(v==1 || v==2 || v==3){
        if(dr_mu<0.2) {
          m_hSvc.h1D("mom_mu_match_ring_vector","","")->Fill(pmomv(v));
          if(prob_ring_angle<0) m_hSvc.h1D("mom_mu_vector_prob_ring_angle_elike_mu_match_ring","","")->Fill(pmomv(v));
          else m_hSvc.h1D("mom_mu_vector_prob_ring_angle_mulike_mu_match_ring","","")->Fill(pmomv(v));
          if(prob_hit<0) m_hSvc.h1D("mom_mu_vector_prob_hit_elike_mu_match_ring","","")->Fill(pmomv(v));
          else m_hSvc.h1D("mom_mu_vector_prob_hit_mulike_mu_match_ring","","")->Fill(pmomv(v));
          if(ip(r)==2) m_hSvc.h1D("mom_mu_vector_pid_elike_mu_match_ring","","")->Fill(pmomv(v));
          else m_hSvc.h1D("mom_mu_vector_pid_mulike_mu_match_ring","","")->Fill(pmomv(v));
          for(int ra=0;ra<n_range_match_mu_mom-1;ra++){
            if(pmomv(v)>range_match_mu_mom[ra] && pmomv(v)<range_match_mu_mom[ra+1]){
              m_hSvc.h1D(Form("prob_ring_angle_mu_match_ringmom%d_%d",range_match_mu_mom[ra],range_match_mu_mom[ra+1]),"","")->Fill(prob_ring_angle);
              m_hSvc.h1D(Form("prob_hit_mu_match_ringmom%d_%d",range_match_mu_mom[ra],range_match_mu_mom[ra+1]),"","")->Fill(prob_hit);
            }
          }
          TVector3 vector_mu = GetTVectorVector(v);
          TVector3 vector_ring = GetTVectorRing(r);
          float dtheta = vector_mu.Theta() - vector_ring.Theta();
          float dphi = vector_mu.DeltaPhi(vector_ring);
          float amomm_residual = (pmomv(v)-amomm(r))/pmomv(v);
          for(int ra=0;ra<n_range_match_mu_mom-1;ra++){
            if(pmomv(v)>range_match_mu_mom[ra] && pmomv(v)<range_match_mu_mom[ra+1]){
              m_hSvc.h1D(Form("dtheta_mu_mom%d_%d",range_match_mu_mom[ra],range_match_mu_mom[ra+1]),"","")->Fill(dtheta);
              m_hSvc.h1D(Form("dphi_mu_mom%d_%d",range_match_mu_mom[ra],range_match_mu_mom[ra+1]),"","")->Fill(dphi);
              m_hSvc.h1D(Form("residual_mulike_mom_mu_mom%d_%d",range_match_mu_mom[ra],range_match_mu_mom[ra+1]),"","")->Fill(amomm_residual);
            }
          }
        }
      }
    }
  }
  m_hSvc.h1D("nMulikeRing","","")->Fill(nRing - n_elike);
  m_hSvc.h1D(Form("nMulikeRing_nring%d",nRing),"","")->Fill(nRing - n_elike);
  m_hSvc.h2D("minMom_maxMom_each_charged_leptons_vector","","")->Fill(min_mom, min_mom_3rd);
  m_hSvc.h1D("minMom_each_charged_leptons_vector","","")->Fill(min_mom);
  m_hSvc.h1D("maxMom_each_charged_leptons_vector","","")->Fill(min_mom_3rd);
  m_hSvc.h1D(Form("minMom_each_charged_leptons_vector_%s",proton_type[is_free_proton].c_str()),"","")->Fill(min_mom);
  m_hSvc.h1D(Form("maxMom_each_charged_leptons_vector_%s",proton_type[is_free_proton].c_str()),"","")->Fill(min_mom_3rd);
  m_hSvc.h1D(Form("minMom_each_charged_leptons_vector_nring%d",nRing),"","")->Fill(min_mom);
  m_hSvc.h1D(Form("middleMom_each_charged_leptons_vector_nring%d",nRing),"","")->Fill(min_mom_2nd);
  m_hSvc.h1D(Form("maxMom_each_charged_leptons_vector_nring%d",nRing),"","")->Fill(min_mom_3rd);

  TLorentzVector vec_mu1,vec_mu2,vec_mu3;
  vec_mu1 = GetTLorentzVectorVector(1);
  vec_mu2 = GetTLorentzVectorVector(2);
  vec_mu3 = GetTLorentzVectorVector(3);
  float mass_mumumu_vector = (vec_mu1+vec_mu2+vec_mu3).M();
  float mom_mumumu_vector = (vec_mu1+vec_mu2+vec_mu3).P();
  m_hSvc.h1D(Form("mass_charged_leptons_vector_%s",proton_type[is_free_proton].c_str()),"","")->Fill(mass_mumumu_vector);
  m_hSvc.h1D(Form("mom_charged_leptons_vector_%s",proton_type[is_free_proton].c_str()),"","")->Fill(mom_mumumu_vector);

  m_hSvc.h1D("nRing","","")->Fill(nring(0));

  //apply selection here

  if(wallv(0)<200) return;

  pass_cut[0][0]=true;
  if ( evis(0) > 30.0 && wall(0) > 200.0 && nhitac(0) <  nhitac_cut[skgen] ){//FC & FV
    pass_cut[1][0]=true;
  }

  if(nring(0)==2) pass_cut[2][0]=true;
  if(nring(0)==3) pass_cut[2][1]=true;
  if(nring(0)==4) pass_cut[2][2]=true;

  if(n_mulike==0) pass_cut[3][0]=true;
  if(n_mulike==1) pass_cut[3][1]=true;
  if(n_mulike==2) pass_cut[3][2]=true;
  if(n_mulike==3) pass_cut[3][3]=true;

  // michel (decay) electron cut
  if(ndcy(0)==0) pass_cut[4][0]=true;
  if(ndcy(0)==1) pass_cut[4][1]=true;
  if(ndcy(0)==2) pass_cut[4][2]=true;
  if(ndcy(0)==3) pass_cut[4][3]=true;

  vector<TLorentzVector> mu_cand;
  mu_cand.clear();
  for(int r=0;r<nRing;r++){
    if(ip(r)==3) mu_cand.push_back(GetTLorentzVectorRing(r,2));
  }

  pass_cut[5][0]=true;//no pi0 mass cut

  TLorentzVector total_vec;
  if(n_mulike==3) total_vec = mu_cand[0] + mu_cand[1] + mu_cand[2];
  if(n_mulike==2) total_vec = mu_cand[0] + mu_cand[1];
  total_mass = total_vec.M();
  total_mom = total_vec.P();
  if(total_mass>800 && total_mass<1050 && total_mom<100){//total mass & low momentum
    pass_cut[6][0]=true;
  }
  if(total_mass>800 && total_mass<1050 && total_mom>100 && total_mom<250){//total mass & high momentum
    pass_cut[7][0]=true;
  }

  TLorentzVector all_ring_vec,all_mulike_vec;
  for(int r=0;r<nRing;r++){
    TLorentzVector temp_vec = 
      (ip(r)==2)? GetTLorentzVectorRing(r,0) : GetTLorentzVectorRing(r,2);
    all_ring_vec += temp_vec;
    all_mulike_vec += GetTLorentzVectorRing(r,2);
  }
  all_ring_mass = all_ring_vec.M();
  all_ring_mom = all_ring_vec.P();
  all_mulike_mass = all_mulike_vec.M();
  all_mulike_mom = all_mulike_vec.P();
  if(kDebugMode) std::cout << "all ring mass/mom=" << all_ring_mass << "/" << all_ring_mom << std::endl;
  if(kDebugMode) std::cout << "all mulike mass/mom=" << all_mulike_mass << "/" << all_mulike_mom << std::endl;

  if(all_mulike_mass>800 && all_mulike_mass<1050 && all_mulike_mom<100){//all_mulike mass & low momentum
    pass_cut[8][0]=true;
  }
  if(all_mulike_mass>800 && all_mulike_mass<1050 && all_mulike_mom>100 && all_mulike_mom<250){//all_mulike mass & high momentum
    pass_cut[9][0]=true;
  }

  MakeCutFlow();

}

void OscNtupleManager::Process_pemumu(){

  for(int v=0;v<nPar;v++){
    int pidgen = ipv(v);
    if(v==1 || v==2 || v==3) {
      m_hSvc.h1D("mom_each_charged_leptons_vector","","")->Fill(pmomv(v));
      TLorentzVector temp_lep = GetTLorentzVectorVector(v);
      m_hSvc.h1D("theta_each_charged_leptons_vector","","")->Fill(temp_lep.Theta());
      m_hSvc.h1D("phi_each_charged_leptons_vector","","")->Fill(temp_lep.Phi());
    }
  }
  TLorentzVector vec_e1,vec_mu1,vec_mu2;
  vec_e1 = GetTLorentzVectorVector(1);
  vec_mu1 = GetTLorentzVectorVector(2);
  vec_mu2 = GetTLorentzVectorVector(3);
  float mass_emumu_vector = (vec_e1+vec_mu1+vec_mu2).M();
  float mom_emumu_vector = (vec_e1+vec_mu1+vec_mu2).P();
  m_hSvc.h1D(Form("mass_charged_leptons_vector_%s",proton_type[is_free_proton].c_str()),"","")->Fill(mass_emumu_vector);
  m_hSvc.h1D(Form("mom_charged_leptons_vector_%s",proton_type[is_free_proton].c_str()),"","")->Fill(mom_emumu_vector);

  //apply selection here

  if(wallv(0)<200) return;

  pass_cut[0][0]=true;
  if ( evis(0) > 30.0 && wall(0) > 200.0 && nhitac(0) <  nhitac_cut[skgen] ){//FC & FV
    pass_cut[1][0]=true;
  }

  //# of cherenkov ring
  if(nring(0)==2) pass_cut[2][0]=true;
  if(nring(0)==3) pass_cut[2][1]=true;
  if(nring(0)==4) pass_cut[2][2]=true;

  //#of mu-like ring
  if(n_mulike==0) pass_cut[3][0]=true;
  if(n_mulike==1) pass_cut[3][1]=true;
  if(n_mulike==2) pass_cut[3][2]=true;
  if(n_mulike==3) pass_cut[3][3]=true;

  // michel (decay) electron cut
  if(ndcy(0)==0) pass_cut[4][0]=true;
  if(ndcy(0)==1) pass_cut[4][1]=true;
  if(ndcy(0)==2) pass_cut[4][2]=true;
  if(ndcy(0)==3) pass_cut[4][3]=true;


  vector<TLorentzVector> mu_cand;
  TLorentzVector e_cand;
  mu_cand.clear();
  for(int r=0;r<nRing;r++){
    if(ip(r)==2) e_cand = GetTLorentzVectorRing(r,1);
    if(ip(r)==3) mu_cand.push_back(GetTLorentzVectorRing(r,2));
  }

  pass_cut[5][0]=true;//no pi0 mass cut

  TLorentzVector total_vec;
  if(n_mulike==2) total_vec = mu_cand[0] + mu_cand[1] + e_cand;
  if(n_mulike==1) total_vec = mu_cand[0] + e_cand;
  total_mass = total_vec.M();
  total_mom = total_vec.P();
  if(total_mass>800 && total_mass<1050 && total_mom<100){//total mass & low momentum
    pass_cut[6][0]=true;
  }
  if(total_mass>800 && total_mass<1050 && total_mom>100 && total_mom<250){//total mass & high momentum
    pass_cut[7][0]=true;
  }

  TLorentzVector all_ring_vec;
  for(int r=0;r<nRing;r++){
    TLorentzVector temp_vec = 
      (ip(r)==2)? GetTLorentzVectorRing(r,0) : GetTLorentzVectorRing(r,2);
    all_ring_vec += temp_vec;
  }
  all_ring_mass = all_ring_vec.M();
  all_ring_mom = all_ring_vec.P();
  if(kDebugMode) std::cout << "all ring mass/mom=" << all_ring_mass << "/" << all_ring_mom << std::endl;

  MakeCutFlow();

}

void OscNtupleManager::Process_pmuee(){

  vector<TLorentzVector> vec_e;
  vec_e.clear();
  TLorentzVector vec_mu;
  for(int v=0;v<nPar;v++){
    int pidgen = ipv(v);
    if(pidgen==2 || pidgen==3) vec_e.push_back(GetTLorentzVectorVector(v));
    if(pidgen==5) vec_mu = GetTLorentzVectorVector(v);
    if(v==1 || v==2 || v==3) {
      m_hSvc.h1D("mom_each_charged_leptons_vector","","")->Fill(pmomv(v));
      TLorentzVector temp_lep = GetTLorentzVectorVector(v);
      m_hSvc.h1D("theta_each_charged_leptons_vector","","")->Fill(temp_lep.Theta());
      m_hSvc.h1D("phi_each_charged_leptons_vector","","")->Fill(temp_lep.Phi());
    }
  }
  TLorentzVector vec_mu1,vec_e1,vec_e2;
  vec_mu1 = GetTLorentzVectorVector(1);
  vec_e1 = GetTLorentzVectorVector(2);
  vec_e2 = GetTLorentzVectorVector(3);
  float mass_muee_vector = (vec_mu1+vec_e1+vec_e2).M();
  float mom_muee_vector = (vec_mu1+vec_e1+vec_e2).P();
  m_hSvc.h1D(Form("mass_charged_leptons_vector_%s",proton_type[is_free_proton].c_str()),"","")->Fill(mass_muee_vector);
  m_hSvc.h1D(Form("mom_charged_leptons_vector_%s",proton_type[is_free_proton].c_str()),"","")->Fill(mom_muee_vector);

  //apply selection here

  if(wallv(0)<200) return;

  pass_cut[0][0]=true;
  if ( evis(0) > 30.0 && wall(0) > 200.0 && nhitac(0) <  nhitac_cut[skgen] ){//FC & FV
    pass_cut[1][0]=true;
  }

  //# of cherenkov ring
  if(nring(0)==2) pass_cut[2][0]=true;
  if(nring(0)==3) pass_cut[2][1]=true;
  if(nring(0)==4) pass_cut[2][2]=true;

  //#of mu-like ring
  if(n_mulike==0) pass_cut[3][0]=true;
  if(n_mulike==1) pass_cut[3][1]=true;
  if(n_mulike==2) pass_cut[3][2]=true;
  if(n_mulike==3) pass_cut[3][3]=true;

  // michel (decay) electron cut
  if(ndcy(0)==0) pass_cut[4][0]=true;
  if(ndcy(0)==1) pass_cut[4][1]=true;
  if(ndcy(0)==2) pass_cut[4][2]=true;
  if(ndcy(0)==3) pass_cut[4][3]=true;

  vector<TLorentzVector> e_cand;
  TLorentzVector mu_cand;
  e_cand.clear();
  for(int r=0;r<nRing;r++){
    if(ip(r)==2) e_cand.push_back(GetTLorentzVectorRing(r,1));
    if(ip(r)==3) mu_cand = GetTLorentzVectorRing(r,2);
  }

  pass_cut[5][0]=true;//no pi0 mass cut

  TLorentzVector total_vec;
  if(n_elike==2) total_vec = e_cand[0] + e_cand[1] + mu_cand;
  if(n_elike==1) total_vec = e_cand[0] + mu_cand;
  total_mass = total_vec.M();
  total_mom = total_vec.P();
  if(total_mass>800 && total_mass<1050 && total_mom<100){//total mass & low momentum
    pass_cut[6][0]=true;
  }
  if(total_mass>800 && total_mass<1050 && total_mom>100 && total_mom<250){//total mass & high momentum
    pass_cut[7][0]=true;
  }

  TLorentzVector all_ring_vec;
  for(int r=0;r<nRing;r++){
    TLorentzVector temp_vec = 
      (ip(r)==2)? GetTLorentzVectorRing(r,0) : GetTLorentzVectorRing(r,2);
    all_ring_vec += temp_vec;
  }
  all_ring_mass = all_ring_vec.M();
  all_ring_mom = all_ring_vec.P();
  if(kDebugMode) std::cout << "all ring mass/mom=" << all_ring_mass << "/" << all_ring_mom << std::endl;

  MakeCutFlow();

}

void OscNtupleManager::ZeroStructure()
{

  os->Zero();
  nDecayE = 0;
  event_type=-1;
  nRing = nring(0);
  nPar = npar(0);
  nPar2 = npar2(0);
  is_free_proton = (pmomv(0)>0.01)? 0 : 1;
  closest_mass_pi0_reco=0.;
  total_mass=0.;
  total_mom=0.;
  all_ring_mass=0.;
  all_ring_mom=0.;

  n_elike=0;
  for(int r=0;r<nRing;r++)
    if(ip(r)==2) n_elike++;

  n_mulike=0;
  for(int r=0;r<nRing;r++) 
    if(ip(r)==3) n_mulike++;

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
  if(kDebugMode) std::cout << "OscNtupleManager::Oscillate" << std::endl;

  return Neighbor3D();

}

double OscNtupleManager::Neighbor3D(){
  if(kDebugMode) std::cout << "OscNtupleManager::Neighbor3D" << std::endl;

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
  //std::cout << "mode=" << mode(0) << std::endl;
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
  //std::cout << "OscNtupleManager::MakeOscillationPlot" << std::endl;
  //std::cout << "interacion_type=" << interaction_type << std::endl;
  if(interaction_type!=-1)  {
    m_hSvc.h1D(Form("nRing_type%d",interaction_type),"","")->Fill(nring(0),live_time_weight*mc_weight);
    m_hSvc.h1D(Form("nRing_type%d_osc",interaction_type),"","")->Fill(nring(0),live_time_weight*osc_weight);
  }
  if(interaction_type>=1 && interaction_type<=10 && interaction_type!=7){//single ring FC event
    //std::cout << "single ring event?" << std::endl;
    //std::cout << "nring=" << nring(0) << std::endl;
    float zenith_angle = -1.*dir(0,2);
    float momentum = (ip(0)==2)? amome(0) : amomm(0);
    m_hSvc.h1D(Form("zenith_angle_type%d",interaction_type),"","")->Fill(zenith_angle,live_time_weight*mc_weight);
    m_hSvc.h1D(Form("zenith_angle_type%d_osc",interaction_type),"","")->Fill(zenith_angle,live_time_weight*osc_weight);
    m_hSvc.h1D(Form("momentum_log10_type%d",interaction_type),"","")->Fill(log10(momentum),live_time_weight*mc_weight);
    m_hSvc.h1D(Form("momentum_log10_type%d_osc",interaction_type),"","")->Fill(log10(momentum),live_time_weight*osc_weight);
    for(int m=0;m<n_range_momentum-1;m++){
      if(momentum>range_momentum[m] && momentum<range_momentum[m+1]){
        m_hSvc.h1D(Form("zenith_angle_type%d_mom%d_%d",interaction_type,range_momentum[m],range_momentum[m+1]),"","")->Fill(zenith_angle,live_time_weight*mc_weight);
        m_hSvc.h1D(Form("zenith_angle_type%d_mom%d_%d_osc",interaction_type,range_momentum[m],range_momentum[m+1]),"","")->Fill(zenith_angle,live_time_weight*osc_weight);
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
        m_hSvc.h1D(Form("zenith_angle_type%d_mom%d_%d",interaction_type,range_momentum[m],range_momentum[m+1]),"","")->Fill(zenith_angle,live_time_weight*mc_weight);
        m_hSvc.h1D(Form("zenith_angle_type%d_mom%d_%d_osc",interaction_type,range_momentum[m],range_momentum[m+1]),"","")->Fill(zenith_angle,live_time_weight*osc_weight);
      }
    }
    //std::cout << "multi ring event?" << std::endl;
    //std::cout << "nring=" << nring(0) << std::endl;

  }
}

void OscNtupleManager::MakeCutFlow(){

  for(int r=0;r<3;r++){// # of ring
    for(int mu=0;mu<4;mu++){// # of mu-like ring
      for(int m=0;m<4;m++){//# of michel electron
        for(int c=0;c<10;c++){//cut
          //std::cout << "r/mu/m/c=" << r << "/" << mu << "/" << m << "/" << c << std::endl;
          //if(c!=2 && c!=3 && c!=4) std::cout << "pass_cut=" << pass_cut[c][0] << std::endl;
          //else if(c==2) std::cout << "pass_cut=" << pass_cut[c][r] << std::endl;
          //else if(c==3) std::cout << "pass_cut=" << pass_cut[c][mu] << std::endl;
          //else if(c==4) std::cout << "pass_cut=" << pass_cut[c][m] << std::endl;
          if(c<6){
            if(!pass_cut[c][0] && c!=2 && c!=3 && c!=4)  break;
            if((!pass_cut[c][r] && c==2) || (!pass_cut[c][mu] && c==3) || (!pass_cut[c][m] && c==4) ) break;
          }
          else if(!pass_cut[c][0]) continue; 
          MakeBasicPlot(c,r,mu,m);
        }
      }
    }
  }

}

void OscNtupleManager::MakeCutFlowValidate(){

  for(int c=0;c<10;c++){//cut
    if(c<6 && !pass_cut[c][0])  break;
    else if(!pass_cut[c][0]) continue; 
    MakeBasicPlot(c,0,0,0);
  }

}

void OscNtupleManager::MakeBasicPlot(int c, int r, int mu, int p){//cut #, michel e cut, 


  if(process_input!="fcdt" || total_mass<800 || total_mass>1050 || total_mom>250){//for blind analysis
    //m_hSvc.h1D("cut_flow","","")->Fill(c,live_time_weight);
    //m_hSvc.h1D("cut_flow_weight_osc","","")->Fill(c,live_time_weight*osc_weight);
    m_hSvc.h1D(Form("cut_flow_nring%d_mulike%d_michel%d",r,mu,p),"","")->Fill(c,live_time_weight*osc_weight);
  }
  //only live time weight
  m_hSvc.h1D(Form("distance_to_wall_cut%d_nring%d_mulike%d_michel%d",c,r,mu,p),"","")->Fill(wall(0),live_time_weight);
  m_hSvc.h1D(Form("visible_energy_cut%d_nring%d_mulike%d_michel%d",c,r,mu,p),"","")->Fill(evis(0),live_time_weight);
  m_hSvc.h1D(Form("nhit_OD_cut%d_nring%d_mulike%d_michel%d",c,r,mu,p),"","")->Fill(nhitac(0),live_time_weight);
  m_hSvc.h1D(Form("nRing_cut%d_nring%d_mulike%d_michel%d",c,r,mu,p),"","")->Fill(nRing,live_time_weight);
  m_hSvc.h1D(Form("nElikeRing_nring%d_cut%d_nring%d_mulike%d_michel%d",nRing,c,r,mu,p),"","")->Fill(n_elike,live_time_weight);
  m_hSvc.h1D(Form("nMulikeRing_nring%d_cut%d_nring%d_mulike%d_michel%d",nRing,c,r,mu,p),"","")->Fill(n_mulike,live_time_weight);
  m_hSvc.h1D(Form("n_michel_electron_cut%d_nring%d_mulike%d_michel%d",c,r,mu,p),"","")->Fill(ndcy(0),live_time_weight);
  //apply oscillation
  m_hSvc.h1D(Form("distance_to_wall_weight_osc_cut%d_nring%d_mulike%d_michel%d",c,r,mu,p),"","")->Fill(wall(0),live_time_weight*osc_weight);
  m_hSvc.h1D(Form("visible_energy_weight_osc_cut%d_nring%d_mulike%d_michel%d",c,r,mu,p),"","")->Fill(evis(0),live_time_weight*osc_weight);
  m_hSvc.h1D(Form("nhit_OD_weight_osc_cut%d_nring%d_mulike%d_michel%d",c,r,mu,p),"","")->Fill(nhitac(0),live_time_weight*osc_weight);
  m_hSvc.h1D(Form("nRing_weight_osc_cut%d_nring%d_mulike%d_michel%d",c,r,mu,p),"","")->Fill(nRing,live_time_weight*osc_weight);
  m_hSvc.h1D(Form("nElikeRing_nring%d_weight_osc_cut%d_nring%d_mulike%d_michel%d",nRing,c,r,mu,p),"","")->Fill(n_elike,live_time_weight*osc_weight);
  m_hSvc.h1D(Form("nMulikeRing_nring%d_weight_osc_cut%d_nring%d_mulike%d_michel%d",nRing,c,r,mu,p),"","")->Fill(n_mulike,live_time_weight*osc_weight);
  m_hSvc.h1D(Form("n_michel_electron_weight_osc_cut%d_nring%d_mulike%d_michel%d",c,r,mu,p),"","")->Fill(ndcy(0),live_time_weight*osc_weight);
  if(n_elike>=2) {
    m_hSvc.h1D(Form("mass_pi0_reco_cut%d_nring%d_mulike%d_michel%d",c,r,mu,p),"","")->Fill(closest_mass_pi0_reco,live_time_weight);
    m_hSvc.h1D(Form("mass_pi0_reco_weight_osc_cut%d_nring%d_mulike%d_michel%d",c,r,mu,p),"","")->Fill(closest_mass_pi0_reco,live_time_weight*osc_weight);
  }
  m_hSvc.h1D(Form("mass_pi0_reco_elike%d_cut%d_nring%d_mulike%d_michel%d",n_elike,c,r,mu,p),"","")->Fill(closest_mass_pi0_reco,live_time_weight);
  m_hSvc.h1D(Form("mass_pi0_reco_elike%d_weight_osc_cut%d_nring%d_mulike%d_michel%d",n_elike,c,r,mu,p),"","")->Fill(closest_mass_pi0_reco,live_time_weight*osc_weight);
  if(process_input!="fcdt" || total_mass<800 || total_mass>1050){//for blind analysis
    m_hSvc.h1D(Form("mass_proton_reco_cut%d_nring%d_mulike%d_michel%d",c,r,mu,p),"","")->Fill(total_mass,live_time_weight);
    m_hSvc.h1D(Form("mass_proton_reco_elike%d_cut%d_nring%d_mulike%d_michel%d",n_elike,c,r,mu,p),"","")->Fill(total_mass,live_time_weight);
    m_hSvc.h1D(Form("mass_proton_reco_mulike%d_cut%d_nring%d_mulike%d_michel%d",n_mulike,c,r,mu,p),"","")->Fill(total_mass,live_time_weight);
    m_hSvc.h1D(Form("mass_proton_reco_weight_mc_cut%d_nring%d_mulike%d_michel%d",c,r,mu,p),"","")->Fill(total_mass,live_time_weight*mc_weight);
    m_hSvc.h1D(Form("mass_proton_reco_weight_osc_cut%d_nring%d_mulike%d_michel%d",c,r,mu,p),"","")->Fill(total_mass,live_time_weight*osc_weight);
    m_hSvc.h1D(Form("mass_proton_reco_nring%d_elike%d_cut%d_nring%d_mulike%d_michel%d",nRing,n_elike,c,r,mu,p),"","")->Fill(total_mass,live_time_weight);
    m_hSvc.h1D(Form("mass_proton_reco_nring%d_mulike%d_cut%d_nring%d_mulike%d_michel%d",nRing,n_mulike,c,r,mu,p),"","")->Fill(total_mass,live_time_weight);
    m_hSvc.h1D(Form("mass_all_ring_reco_nring%d_cut%d_nring%d_mulike%d_michel%d",nRing,c,r,mu,p),"","")->Fill(all_ring_mass,live_time_weight*osc_weight);
    m_hSvc.h1D(Form("mass_all_mulike_reco_nring%d_cut%d_nring%d_mulike%d_michel%d",nRing,c,r,mu,p),"","")->Fill(all_mulike_mass,live_time_weight*osc_weight);
  }
  if(process_input!="fcdt" || total_mom>250){//for blind analysis
    m_hSvc.h1D(Form("mom_proton_reco_cut%d_nring%d_mulike%d_michel%d",c,r,mu,p),"","")->Fill(total_mom,live_time_weight);
    m_hSvc.h1D(Form("mom_proton_reco_elike%d_cut%d_nring%d_mulike%d_michel%d",n_elike,c,r,mu,p),"","")->Fill(total_mom,live_time_weight);
    m_hSvc.h1D(Form("mom_proton_reco_mulike%d_cut%d_nring%d_mulike%d_michel%d",n_mulike,c,r,mu,p),"","")->Fill(total_mom,live_time_weight);
    m_hSvc.h1D(Form("mom_proton_reco_weight_mc_cut%d_nring%d_mulike%d_michel%d",c,r,mu,p),"","")->Fill(total_mom,live_time_weight*mc_weight);
    m_hSvc.h1D(Form("mom_proton_reco_weight_osc_cut%d_nring%d_mulike%d_michel%d",c,r,mu,p),"","")->Fill(total_mom,live_time_weight*osc_weight);
    m_hSvc.h1D(Form("mom_all_ring_reco_nring%d_cut%d_nring%d_mulike%d_michel%d",nRing,c,r,mu,p),"","")->Fill(all_ring_mom,live_time_weight*osc_weight);
    m_hSvc.h1D(Form("mom_all_mulike_reco_nring%d_cut%d_nring%d_mulike%d_michel%d",nRing,c,r,mu,p),"","")->Fill(all_mulike_mom,live_time_weight*osc_weight);
  }
  if(process_input!="fcdt" || total_mass<800 || total_mass>1050 || total_mom>250){//for blind analysis
    m_hSvc.h2D(Form("mass_mom_proton_reco_cut%d_nring%d_mulike%d_michel%d",c,r,mu,p),"","")->Fill(total_mass,total_mom,live_time_weight);
    m_hSvc.h2D(Form("mass_mom_proton_reco_elike%d_cut%d_nring%d_mulike%d_michel%d",n_elike,c,r,mu,p),"","")->Fill(total_mass,total_mom,live_time_weight);
    m_hSvc.h2D(Form("mass_mom_proton_reco_mulike%d_cut%d_nring%d_mulike%d_michel%d",n_mulike,c,r,mu,p),"","")->Fill(total_mass,total_mom,live_time_weight);
    m_hSvc.h2D(Form("mass_mom_proton_reco_weight_mc_cut%d_nring%d_mulike%d_michel%d",c,r,mu,p),"","")->Fill(total_mass,total_mom,live_time_weight*mc_weight);
    m_hSvc.h2D(Form("mass_mom_proton_reco_weight_osc_cut%d_nring%d_mulike%d_michel%d",c,r,mu,p),"","")->Fill(total_mass,total_mom,live_time_weight*osc_weight);
    m_hSvc.h2D(Form("all_ring_mass_mom_proton_reco_cut%d_nring%d_mulike%d_michel%d",c,r,mu,p),"","")->Fill(all_ring_mass,all_ring_mom,live_time_weight*osc_weight);
    m_hSvc.h2D(Form("all_mulike_mass_mom_proton_reco_cut%d_nring%d_mulike%d_michel%d",c,r,mu,p),"","")->Fill(all_mulike_mass,all_mulike_mom,live_time_weight*osc_weight);
  }

}

int OscNtupleManager::GetParType(int ring_id){
  float min_dr=99.;
  int min_pid=0;
  for(int p=1;p<npar(0);p++){//initial proton is skipped
    float dr = CalcDrRingVector(ring_id,p);
    cout << "pid=" << ipv(p) << endl;
    cout << "dr=" << dr << endl;
    if(/*dr<0.2 &&*/ dr<min_dr){
      min_dr = dr;
      min_pid = ipv(p);
    }
  }
  return min_pid;

}

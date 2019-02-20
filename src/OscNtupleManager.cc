#include "src/OscNtupleManager.h"
#include <fstream>

#include "TLorentzVector.h"
#include "TRandom.h"

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

  ehit_cut_1[0] =  60; 
  ehit_cut_1[1] =  30; 
  ehit_cut_1[2] =  60; 
  ehit_cut_1[3] =  60; 

  ehit_cut_2[0] =  40; 
  ehit_cut_2[1] =  20; 
  ehit_cut_2[2] =  40; 
  ehit_cut_2[3] =  40; 

  llMGMRE       = -1e7 ;
  llMGMRENUE    = -1e7 ;
  llMGMRENUEBar = -1e7 ;

  kUseFiTQun    = false ;
  kUseTauNN     = false ; 
  kDebugMode     = false ; 
  kAllHist     = false ; 
  kMakeNtuple     = false ; 
  kCorrelatedDecay = -1;
  kFermiMotion = -1;
  kOutsideSR = -1;
  kSystNtag = -1;
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
  sample_num=0;event_num=-1;

  for(int r=0;r<r_max;r++){
    for(int mu=0;mu<mu_max;mu++){
      for(int m=0;m<m_max;m++){
        total_with_ntag[r][mu][m]=0;
        lowerSR_without_ntag[r][mu][m]=0;
        higherSR_without_ntag[r][mu][m]=0;
      }
    }
  }

}


OscNtupleManager::~OscNtupleManager( )
{
  delete os;
  delete llm;


  // delete llmgmre;
  // delete llpi0;
}

void OscNtupleManager::CreateBranch(){
  if(kDebugMode) cout << "CreateBranch" << endl;

  //ofile = new TFile("aho.root","recreate");
  ofile = new TFile(output_ntuple.c_str(),"recreate");
  otree = new TTree("osc_tuple", "tree build for SK osc analyses"); 
  otree->Branch("ipnu"   , &o_ipnu    , "ipnu/I"    );
  otree->Branch("itype"   , &o_itype    , "itype/I"    );
  otree->Branch("pnu"   , &o_pnu    , "pnu/F"    );
  otree->Branch("dir"    , &o_dir     , "dir[3]/F"  );
  otree->Branch("amom"   , &o_amom    , "amom/F"    );
  otree->Branch("mode"   , &o_mode    , "mode/I"    );
  otree->Branch("nring"   , &o_nring    , "nring/I"    );
  otree->Branch("nmulike"   , &o_nmulike    , "nmulike/I"    );
  otree->Branch("total_mass"   , &o_total_mass    , "total_mass/F"    );
  otree->Branch("total_mom"   , &o_total_mom    , "total_mom/F"    );
  otree->Branch("dlfct"   , &o_dlfct    , "dlfct/F"    );
  otree->Branch("prob_angle"    , &o_prob_angle     , "prob_angle[5]/F"  );
  otree->Branch("probms_e"    , &o_probms_e     , "probms_e[5]/F"  );
  otree->Branch("probms_mu"    , &o_probms_mu     , "probms_mu[5]/F"  );
  otree->Branch("prob_pattern"    , &o_prob_pattern     , "prob_pattern[5]/F"  );
  otree->Branch("prmslg_e"    , &o_prmslg_e     , "prmslg_e[5]/F"  );
  otree->Branch("prmslg_mu"    , &o_prmslg_mu     , "prmslg_mu[5]/F"  );
  otree->Branch("mmom"    , &o_mmom     , "mmom[5]/F"  );
  otree->Branch("dir_x"   , &o_dir_x    , "dir_x[5]/F"    );
  otree->Branch("dir_y"   , &o_dir_y    , "dir_y[5]/F"    );
  otree->Branch("dir_z"   , &o_dir_z    , "dir_z[5]/F"    );
  otree->Branch("ang"    , &o_ang     , "ang[5]/F"  );
  otree->Branch("ange"    , &o_ange     , "ange[5]/F"  );
  otree->Branch("angm"    , &o_angm     , "angm[5]/F"  );
  otree->Branch("mmom_min"    , &o_mmom_min     , "mmom_min/F"  );
  otree->Branch("mmom_mid"    , &o_mmom_mid     , "mmom_mid/F"  );
  otree->Branch("mmom_max"    , &o_mmom_max     , "mmom_max/F"  );
  otree->Branch("vertex_x"   , &o_vertex_x    , "vertex_x/F"    );
  otree->Branch("vertex_y"   , &o_vertex_y    , "vertex_y/F"    );
  otree->Branch("vertex_z"   , &o_vertex_z    , "vertex_z/F"    );

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
    m_hSvc.create1D("hstate","",5,0,5);
    m_hSvc.create1D("fermi_momentum","",15,0,300);
    m_hSvc.create1D("fermi_dirx","",20,-1,1);
    m_hSvc.create1D("fermi_diry","",20,-1,1);
    m_hSvc.create1D("fermi_dirz","",20,-1,1);
    m_hSvc.create1D("nRing","",10,0,10);
    m_hSvc.create1D("nMulikeRing_angle","",6,0,6);
    m_hSvc.create1D("ring_counting_likelihood","",100,-15,15);
    m_hSvc.create1D("true_dirx_lepton","",100,-1,1);
    m_hSvc.create1D("true_diry_lepton","",100,-1,1);
    m_hSvc.create1D("true_dirz_lepton","",100,-1,1);
    m_hSvc.create1D("true_mom_electron","",40,0,800);
    m_hSvc.create1D("true_mom_muon","",40,0,800);
    m_hSvc.create1D("diff_vertex_x","",100,-100,100);
    m_hSvc.create1D("diff_vertex_y","",100,-100,100);
    m_hSvc.create1D("diff_vertex_z","",100,-100,100);
    m_hSvc.create1D("diff_vertex_r","",100,-100,100);
    m_hSvc.create1D("true_angle_lepton_and_ring","",36,0,180);
    m_hSvc.create1D("true_mom_electron_match_ring","",40,0,800);
    m_hSvc.create1D("true_mom_muon_match_ring","",40,0,800);
    m_hSvc.create1D("true_mom_electron_match_ring_angle_elike","",40,0,800);
    m_hSvc.create1D("true_mom_muon_match_ring_angle_mulike","",40,0,800);
    m_hSvc.create1D("true_angle_min_mid_lepton","",36,0,180);
    m_hSvc.create1D("true_angle_min_max_lepton","",36,0,180);
    m_hSvc.create1D("true_angle_mid_max_lepton","",36,0,180);
    for(int r=2;r<4;r++){
      m_hSvc.create1D(Form("true_angle_lepton_and_ring_nring%d",r),"",36,0,180);
      m_hSvc.create1D(Form("true_mom_lepton_nring%d",r),"",40,0,800);
      m_hSvc.create1D(Form("true_mom_muon_nring%d",r),"",40,0,800);
      m_hSvc.create1D(Form("true_mom_muon_match_ring_nring%d",r),"",40,0,800);
      m_hSvc.create1D(Form("true_mom_muon_match_ring_angle_mulike_nring%d",r),"",40,0,800);
      m_hSvc.create1D(Form("true_mom_electron_match_ring_nring%d",r),"",40,0,800);
      m_hSvc.create1D(Form("true_mom_electron_match_ring_angle_elike_nring%d",r),"",40,0,800);
      for(int m1=0;m1<4;m1++){
        m_hSvc.create1D(Form("true_mom_muon_nring%d_mulike%d",r,m1),"",40,0,800);
        for(int f=0;f<2;f++){
          m_hSvc.create1D(Form("total_mass_nring%d_mulike%d_fp%d",r,m1,f),"",50,750,1050);
          m_hSvc.create1D(Form("residual_total_mass_nring%d_mulike%d_fp%d",r,m1,f),"",80,-0.2,0.2);
          m_hSvc.create1D(Form("diff_total_mass_nring%d_mulike%d_fp%d",r,m1,f),"",100,-100,100);
          m_hSvc.create1D(Form("total_mom_nring%d_mulike%d_fp%d",r,m1,f),"",100,0,1000);
          m_hSvc.create1D(Form("residual_total_mom_nring%d_mulike%d_fp%d",r,m1,f),"",80,-0.2,0.2);
          m_hSvc.create1D(Form("diff_total_mom_nring%d_mulike%d_fp%d",r,m1,f),"",100,-100,100);
        }
      }
      for(int p=0;p<6;p++){
        m_hSvc.create1D(Form("residual_emom_mom%d_%d_nring%d",100*p,100+100*p,r),"",80,-0.4,0.4);
        m_hSvc.create1D(Form("residual_mmom_mom%d_%d_nring%d",100*p,100+100*p,r),"",80,-0.4,0.4);
      }
    }
    m_hSvc.create1D("prob_angle_electron","",100,-10,10);
    m_hSvc.create1D("prob_angle_muon","",100,-10,10);
    for(int p=0;p<6;p++){
      m_hSvc.create1D(Form("prob_angle_electron_mom%d_%d",100*p,100+100*p),"",100,-10,10);
      m_hSvc.create1D(Form("prob_angle_muon_mom%d_%d",100*p,100+100*p),"",100,-10,10);
      m_hSvc.create1D(Form("residual_emom_mom%d_%d",100*p,100+100*p),"",80,-0.4,0.4);
      m_hSvc.create1D(Form("residual_mmom_mom%d_%d",100*p,100+100*p),"",80,-0.4,0.4);
    }

    for(int f=0;f<2;f++){
      m_hSvc.create1D(Form("true_mom_lepton_fp%d",is_free_proton),"",40,0,800);
      m_hSvc.create1D(Form("total_true_mass_fp%d",f),"",125,0,1250);
      m_hSvc.create1D(Form("total_true_mom_fp%d",f),"",100,0,1000);
      m_hSvc.create1D(Form("true_min_mom_lepton_fp%d",f),"",40,0,800);
      m_hSvc.create1D(Form("true_mid_mom_lepton_fp%d",f),"",40,0,800);
      m_hSvc.create1D(Form("true_max_mom_lepton_fp%d",f),"",40,0,800);
      m_hSvc.create1D(Form("angle_min_mom_lepton_proton_fp%d",f),"",36,0,180);
      m_hSvc.create1D(Form("angle_mid_mom_lepton_proton_fp%d",f),"",36,0,180);
      m_hSvc.create1D(Form("angle_max_mom_lepton_proton_fp%d",f),"",36,0,180);
    }
  }
  if(process_mode=="subgev_multiring"){
    int dwall_thr[3] = {50,100,150};
    for(int d=0;d<3;d++){
      m_hSvc.create1D(Form("distance_to_wall_thr%d",dwall_thr[d]),"",36,0,1800);
      for(int r=2;r<6;r++){
        m_hSvc.create1D(Form("distance_to_wall_thr%d_nring%d",dwall_thr[d],r),"",36,0,1800);
      }
    }
    for(int m=1;m<4;m++) m_hSvc.create1D(Form("n_michel_electron_mulike%d",m),"",5,0,5);

    //check background component
    for (int m=1;m<55;m++){
      for(int t=0;t<2;t++){
        for(int mu=1;mu<4;mu++){
          m_hSvc.create1D(Form("n_michel_electron_mulike%d_type%d_mode_pos%d",mu,t,m),"",5,0,5);
          m_hSvc.create1D(Form("n_michel_electron_mulike%d_type%d_mode_neg%d",mu,t,m),"",5,0,5);
        }
      }
    }
  }

  if(process_mode=="subgev_onemulike" || process_mode=="subgev_oneelike"){
   
    m_hSvc.create1D("n_michel_electron","",5,0,5);
    m_hSvc.create1D("nRing","",5,0,5);
    m_hSvc.create1D("nMulikeRing_angle","",5,0,5);

    //check background component
    for (int m=1;m<55;m++){
      for(int t=0;t<2;t++){
        m_hSvc.create1D(Form("n_michel_electron_type%d_mode_pos%d",t,m),"",5,0,5);
        m_hSvc.create1D(Form("n_michel_electron_type%d_mode_neg%d",t,m),"",5,0,5);
      }
    }

  }

  r_max=3;
  int cut_max=10;
  if(process_mode.find("p_epi")!=std::string::npos || process_mode.find("p_eee")!=std::string::npos){
    mu_max=1;m_max=1;}
  if(process_mode.find("p_mupi")!=std::string::npos || process_mode.find("p_muee")!=std::string::npos
      || process_mode.find("p_eemu")!=std::string::npos){
    mu_max=1;m_max=3;}
  if(process_mode.find("p_mumumu")!=std::string::npos){
    mu_max=1;m_max=4;}
  if(process_mode.find("p_emumu")!=std::string::npos || process_mode.find("p_mumue")!=std::string::npos){
    mu_max=1;m_max=3;}

  if(kAllHist && !kCheckBkg){
    for(int c=0;c<cut_max;c++){//cut type
      for(int r=0;r<r_max;r++){//cut pattern of nring
        for(int mu=0;mu<mu_max;mu++){//cut pattern of mu-like ring
          for(int p=0;p<m_max;p++){//cut pattern of michel electron
            if(kSystNtag){
              m_hSvc.create1D(Form("ntag_multiplicity_cut%d_nring%d_mulike%d_michel%d",c,r,mu,p),"",10,0,10);
              for(int e=0;e<11;e++){
                int eff = e*10;
                m_hSvc.create1D(Form("n_tagged_neutron_exp_eff%d_cut%d_nring%d_mulike%d_michel%d",eff,c,r,mu,p),"",10,0,10);
              }
              continue;
            }
            m_hSvc.create1D(Form("n_true_neutron_cut%d_nring%d_mulike%d_michel%d",c,r,mu,p),"",15,0,15);
            m_hSvc.create1D(Form("interaction_mode_cut%d_nring%d_mulike%d_michel%d",c,r,mu,p),"",92,1,93);
            m_hSvc.create1D(Form("distance_to_wall_cut%d_nring%d_mulike%d_michel%d",c,r,mu,p),"",36,0,1800);
            m_hSvc.create1D(Form("visible_energy_cut%d_nring%d_mulike%d_michel%d",c,r,mu,p),"",100,0,500);
            m_hSvc.create1D(Form("nhit_OD_cut%d_nring%d_mulike%d_michel%d",c,r,mu,p),"",100,0,100);
            m_hSvc.create1D(Form("nRing_cut%d_nring%d_mulike%d_michel%d",c,r,mu,p),"",10,0,10);
            m_hSvc.create1D(Form("n_michel_electron_cut%d_nring%d_mulike%d_michel%d",c,r,mu,p),"",5,0,5);
            m_hSvc.create1D(Form("ntag_multiplicity_cut%d_nring%d_mulike%d_michel%d",c,r,mu,p),"",10,0,10);
            m_hSvc.create1D(Form("mass_pi0_reco_cut%d_nring%d_mulike%d_michel%d",c,r,mu,p),"",125,0,1250);
            m_hSvc.create1D(Form("mass_proton_reco_cut%d_nring%d_mulike%d_michel%d",c,r,mu,p),"",125,0,1250);
            m_hSvc.create1D(Form("mom_proton_reco_cut%d_nring%d_mulike%d_michel%d",c,r,mu,p),"",100,0,1000);
            m_hSvc.create1D(Form("mass_two_elike_reco_cut%d_nring%d_mulike%d_michel%d",c,r,mu,p),"",125,0,1250);
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
              m_hSvc.create1D(Form("distance_to_wall_cut%d_nring%d_mulike%d_michel%d_fp%d",c,r,mu,p,f),"",36,0,1800);
              m_hSvc.create1D(Form("visible_energy_cut%d_nring%d_mulike%d_michel%d_fp%d",c,r,mu,p,f),"",100,0,500);
              m_hSvc.create1D(Form("nhit_OD_cut%d_nring%d_mulike%d_michel%d_fp%d",c,r,mu,p,f),"",100,0,100);
              m_hSvc.create1D(Form("nRing_cut%d_nring%d_mulike%d_michel%d_fp%d",c,r,mu,p,f),"",10,0,10);
              m_hSvc.create1D(Form("n_michel_electron_cut%d_nring%d_mulike%d_michel%d_fp%d",c,r,mu,p,f),"",5,0,5);
              m_hSvc.create1D(Form("ntag_multiplicity_cut%d_nring%d_mulike%d_michel%d_fp%d",c,r,mu,p,f),"",10,0,10);
              m_hSvc.create1D(Form("mass_pi0_reco_cut%d_nring%d_mulike%d_michel%d_fp%d",c,r,mu,p,f),"",125,0,1250);
              m_hSvc.create1D(Form("mass_proton_reco_cut%d_nring%d_mulike%d_michel%d_fp%d",c,r,mu,p,f),"",125,0,1250);
              m_hSvc.create1D(Form("mom_proton_reco_cut%d_nring%d_mulike%d_michel%d_fp%d",c,r,mu,p,f),"",100,0,1000);
              m_hSvc.create1D(Form("mass_two_elike_reco_cut%d_nring%d_mulike%d_michel%d_fp%d",c,r,mu,p,f),"",125,0,1250);
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
        for(int m=0;m<m_max;m++){
          m_hSvc.create1D(Form("cut_flow_nring%d_mulike%d_michel%d",r,mu,m),"",11,0,11);
          for(int f=0;f<2;f++){
            m_hSvc.create1D(Form("cut_flow_nring%d_mulike%d_michel%d_fp%d",r,mu,m,f),"",11,0,11);
          }
        }
      }
    }

  }//all hist

  if(kCheckBkg){//default is fcmc
    for(int r=0;r<r_max;r++){//cut pattern of nring
      for(int mu=0;mu<mu_max;mu++){//cut pattern of mu-like ring
        for(int p=0;p<m_max;p++){//cut pattern of michel electron
          for (int m=1;m<55;m++){
            for(int t=0;t<2;t++){
              m_hSvc.create1D(Form("cut_flow_nring%d_mulike%d_michel%d_type%d_mode_pos%d",r,mu,p,t,m),"",11,0,11);
              m_hSvc.create1D(Form("cut_flow_nring%d_mulike%d_michel%d_type%d_mode_neg%d",r,mu,p,t,m),"",11,0,11);
              for(int c=0;c<10;c++){//cut type
                m_hSvc.create1D(Form("nRing_cut%d_nring%d_mulike%d_michel%d_type%d_mode_pos%d",c,r,mu,p,t,m),"",10,0,10);
                m_hSvc.create1D(Form("nRing_cut%d_nring%d_mulike%d_michel%d_type%d_mode_neg%d",c,r,mu,p,t,m),"",10,0,10);
                m_hSvc.create1D(Form("nElikeRing_angle_nring2_cut%d_nring%d_mulike%d_michel%d_type%d_mode_pos%d",c,r,mu,p,t,m),"",6,0,6);
                m_hSvc.create1D(Form("nElikeRing_angle_nring3_cut%d_nring%d_mulike%d_michel%d_type%d_mode_pos%d",c,r,mu,p,t,m),"",6,0,6);
                m_hSvc.create1D(Form("nElikeRing_angle_nring2_cut%d_nring%d_mulike%d_michel%d_type%d_mode_neg%d",c,r,mu,p,t,m),"",6,0,6);
                m_hSvc.create1D(Form("nElikeRing_angle_nring3_cut%d_nring%d_mulike%d_michel%d_type%d_mode_neg%d",c,r,mu,p,t,m),"",6,0,6);
                m_hSvc.create1D(Form("nMulikeRing_angle_nring2_cut%d_nring%d_mulike%d_michel%d_type%d_mode_pos%d",c,r,mu,p,t,m),"",6,0,6);
                m_hSvc.create1D(Form("nMulikeRing_angle_nring3_cut%d_nring%d_mulike%d_michel%d_type%d_mode_pos%d",c,r,mu,p,t,m),"",6,0,6);
                m_hSvc.create1D(Form("nMulikeRing_angle_nring2_cut%d_nring%d_mulike%d_michel%d_type%d_mode_neg%d",c,r,mu,p,t,m),"",6,0,6);
                m_hSvc.create1D(Form("nMulikeRing_angle_nring3_cut%d_nring%d_mulike%d_michel%d_type%d_mode_neg%d",c,r,mu,p,t,m),"",6,0,6);
                m_hSvc.create1D(Form("n_michel_electron_cut%d_nring%d_mulike%d_michel%d_type%d_mode_pos%d",c,r,mu,p,t,m),"",5,0,5);
                m_hSvc.create1D(Form("n_michel_electron_cut%d_nring%d_mulike%d_michel%d_type%d_mode_neg%d",c,r,mu,p,t,m),"",5,0,5);
                m_hSvc.create1D(Form("ntag_multiplicity_cut%d_nring%d_mulike%d_michel%d_type%d_mode_pos%d",c,r,mu,p,t,m),"",10,0,10);
                m_hSvc.create1D(Form("ntag_multiplicity_cut%d_nring%d_mulike%d_michel%d_type%d_mode_neg%d",c,r,mu,p,t,m),"",10,0,10);
                m_hSvc.create1D(Form("mass_proton_reco_cut%d_nring%d_mulike%d_michel%d_type%d_mode_pos%d",c,r,mu,p,t,m),"",125,0,1250);
                m_hSvc.create1D(Form("mass_proton_reco_cut%d_nring%d_mulike%d_michel%d_type%d_mode_neg%d",c,r,mu,p,t,m),"",125,0,1250);
                m_hSvc.create1D(Form("mass_two_elike_reco_cut%d_nring%d_mulike%d_michel%d_type%d_mode_pos%d",c,r,mu,p,t,m),"",125,0,1250);
                m_hSvc.create1D(Form("mass_two_elike_reco_cut%d_nring%d_mulike%d_michel%d_type%d_mode_neg%d",c,r,mu,p,t,m),"",125,0,1250);
                m_hSvc.create1D(Form("mom_proton_reco_cut%d_nring%d_mulike%d_michel%d_type%d_mode_pos%d",c,r,mu,p,t,m),"",100,0,1000);
                m_hSvc.create1D(Form("mom_proton_reco_cut%d_nring%d_mulike%d_michel%d_type%d_mode_neg%d",c,r,mu,p,t,m),"",100,0,1000);
              }
            }
          }
        }
      }
    }
  }

}

void OscNtupleManager::LoadWrappers()
{

  dm->Get("hstate"    , hstate     );
  dm->Get( "pnu"  , pnu    );
  dm->Get( "ndcy"  , ndcy    );
  dm->Get( "ngate"  , ngate    );
  dm->Get( "nbye"  , nbye    );
  dm->Get("nmue"    , nmue);
  dm->Get("evis"  , evis   );
  dm->Get("etime"   , etime);
  dm->Get("etype"   , etype);
  dm->Get("ehit"    , ehit);
  dm->Get("egood"   , egood);
  dm->Get("wall"  , wall   );
  dm->Get("nring" , nring  );
  dm->Get("numnu" , numnu  );
  dm->Get("npar" , npar  );
  dm->Get("npar2" , npar2  );
  dm->Get("Npvc" , Npvc  );
  dm->Get("Ipvc" , Ipvc  );
  dm->Get("Ichvc" , Ichvc  );
  dm->Get("Abspvc" , Abspvc  );
  dm->Get("pscnd" , pscnd  );
  dm->Get("nscndprt" , nscndprt  );
  dm->Get("iprtscnd" , iprtscnd  );
  dm->Get("iprntprt" , iprntprt  );
  dm->Get("lmecscnd" , lmecscnd  );
  dm->Get("iprnttrk" , iprnttrk  );
  dm->Get("Dlfct" , Dlfct  );
  dm->Get("Iorgvc" , Iorgvc  );
  dm->Get("Iflvc" , Iflvc  );
  dm->Get("iprntidx" , iprntidx   );
  dm->Get("nhitac", nhitac );
  dm->Get("nev", nev );
  dm->Get("nsub", nsub );
  dm->Get("nsube", nsube );
  dm->Get("date", date );
  dm->Get("time", time );
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
  cout << "kDebugMode=" << kDebugMode << endl;
  cout << "kCorrelatedDecay=" << kCorrelatedDecay << endl;
  cout << "kFermiMotion=" << kFermiMotion << endl;
  cout << "kOutsideSR=" << kOutsideSR << endl;
  cout << "kSystNtag=" << kSystNtag << endl;

  TFile *f_fm;
  TH1 *h_fm=new TH1F();
  if(kFermiMotion){//fermi gas model
    f_fm = new TFile("make_hist/output/ratio_fermi_motion.root");
    h_fm = (TH1F*) f_fm->Get("ratio_fermi_motion");
    cout << "h_fm entries=" << h_fm->GetEntries() << endl;
  }

  for( int i = startEntry ; i < endEntry ; i ++ )
  {
    if(kDebugMode) cout << "Entry is " << i << endl;
    else if(i%10000==0) 
      cout << "Entry is " << i << endl;

    _ltree->GetEntry(i-1);//just for culculating sample number
    int nev_before = nev(0);
    _ltree->GetEntry(i);

    ZeroStructure();

    /*if(kMakeNtuple){
      FillNtuple();
      continue;
    }*/

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
    if(kCorrelatedDecay==0){//no correlated decay
      if(hstate(0)>=1 && hstate(0)<=3) weight = 20./18;//adjust the ratio of free proton   
      if(hstate(0)==4) weight = 0;    
    }
    if(kCorrelatedDecay==2){//20% correlated decay
      if(hstate(0)>=1 && hstate(0)<=3) weight = 16./18;//adjust the ratio of free proton   
      if(hstate(0)==4) weight *= 2;    
    }
    if(kFermiMotion && hstate(0)!=0 && hstate(0)!=4){//fermi gas model
      if(kDebugMode) cout << "proton momentum=" << pmomv(0) << endl;
      for(int p=0;p<h_fm->GetNbinsX();p++){
        float xmin = h_fm->GetBinLowEdge(p+1);
        float xmax = xmin + h_fm->GetBinWidth(p+1);
        if(kDebugMode) cout << "p/xmin/xmax=" << p << "/" << xmin << "/" << xmax << endl;
        if(pmomv(0)>xmin && pmomv(0)<xmax){//range of proton momentum
          if(kDebugMode) cout << "ratio=" << h_fm->GetBinContent(p+1) << endl;
          weight *= h_fm->GetBinContent(p+1);
          break;
        }
      }
    }
    
    if(nev(0)<nev_before) {
      sample_num++;
      event_num=-1;
    }
    event_num++;
    if(kDebugMode) {
      cout << "num sample/evt=" << sample_num << "/" << event_num << endl;
      cout << "nsub/nev=" << nsub(0) << "/" << nev(0) << endl;
      cout << "nsube/ndcy/ngate/nbye=" << nsube(0) << "/" << ndcy(0) << "/" << ngate(0) << "/" << nbye(0) << endl;
      cout << "event date " << date(0) << "/" << date(1) << "/" << date(2) << endl;
      cout << "event time " << time(0) << "/" << time(1) << "/" << time(2) << endl;
      cout << "mode/nutype=" << mode(0) << "/" << ipnu(0) << endl;
      cout << "nring/n_elike/n_mulike=" << nRing << "/" << n_elike_angle << "/" << n_mulike_angle << endl;
      cout << "nDecayE=" << nDecayE << endl;
      cout << "d_wall true/reco=" << wallv(0) << "/" << wall(0) << endl;
      cout << "true vertex position x/y/z=" << posv(0) << "/" << posv(1) << "/" << posv(2) << endl;
      cout << "reco vertex position x/y/z=" << pos(0) << "/" << pos(1) << "/" << pos(2) << endl;
      cout << "live time weight is " << live_time_weight << endl;
      cout << "weight mc is " << mc_weight << endl;
      cout << "weight osc is " << osc_weight << endl;
      cout << "weight=" << weight << endl;
      cout << "npar/npar2=" << npar(0) << "/" << npar2(0) << endl;
      cout << "Dlfct=" << Dlfct(0) << endl;
      if(process_input!="fcmc" && process_input!="fcdt") cout << "hstate=" << hstate(0) << endl;
      if(kCheckBkg) cout << "Npvc/daughter=" << Npvc(0) << "/" << nscndprt(0) << endl;
      cout << "numnu=" << numnu(0) << endl;
      cout << "### vector list at detector simulation ###" << endl;
      for(int m=0;m<npar(0);m++){
        cout << "particle_" << m+1 << " pid=" << ipv(m) << " mom=" << pmomv(m) 
          << " dir x/y/z=" << dirv(m,0) << "/" << dirv(m,1) << "/" << dirv(m,2) << endl;
      }
      cout << "### vector list of secondaries at detector simulation ###" << endl;
      for(int t=0;t<npar2(0);t++){
        cout << "daughter particle from " << iorg(t) << " pid=" << ipv2(t) << " mom=" << pmomv2(t)
          << " dir x/y/z=" << dirv2(t,0) << "/" << dirv2(t,1) << "/" << dirv2(t,2) << endl;
      }
      cout << "### vector list at neutrino interaction ###" << endl;
      for(int n=0;n<numnu(0);n++){
        cout << "numnu_" << n+1 << " ipnu=" << ipnu(n) << endl;
      }
      if(kCheckBkg){
        cout << "### Copy of VECT and NEWORK primary stacks with additional information ###" << endl;
        for(int v=0;v<Npvc(0);v++){
          cout << "pvc_" << v+1 << " pid=" << Ipvc(v) << " org=" << Iorgvc(v) << " mom=" << Abspvc(v) << " flv=" << Iflvc(v) << endl;
        }
        cout << "### Copy of VECT2 secondaries stack with additional information ###" << endl;
        for(int d=0;d<nscndprt(0);d++){
          float momx = pscnd(d,0);
          float momy = pscnd(d,1);
          float momz = pscnd(d,2);
          float mom = sqrt(momx*momx+momy*momy+momz*momz);
          cout << "daughter_" << d+1 << " pid=" << iprtscnd(d) << " org iprntprt/iprnttrk/iprntidx=" << iprntprt(d) << "/" << iprnttrk(d) << "/" << iprntidx(d) 
            << " lmecscnd=" << lmecscnd(d) << " mom=" << mom << endl;
        }
      }
    }

    if(process_mode.find("p_epi")!=std::string::npos) Process_pepi();
    if(process_mode.find("p_mupi")!=std::string::npos) Process_pmupi();
    if(process_mode.find("p_eee")!=std::string::npos) Process_peee();
    if(process_mode.find("p_muee")!=std::string::npos) Process_pmuee();
    if(process_mode.find("p_eemu")!=std::string::npos) Process_pmuee();
    if(process_mode.find("p_emumu")!=std::string::npos) Process_pemumu();
    if(process_mode.find("p_mumue")!=std::string::npos) Process_pemumu();
    if(process_mode.find("p_mumumu")!=std::string::npos) Process_pmumumu();
    if(process_mode.find("single_e")!=std::string::npos) Process_single();
    if(process_mode.find("single_mu")!=std::string::npos) Process_single();
    if(process_mode.find("subgev_multiring")!=std::string::npos) Process_subgev_multiring();
    if(process_mode.find("subgev_onemulike")!=std::string::npos) Process_subgev_onemulike();
    if(process_mode.find("subgev_oneelike")!=std::string::npos) Process_subgev_oneelike();

  }

  cout << " Classification: " << endl
    << "   Written: " << Good << endl
    << "   Skipped: " << Bad << endl
    << endl;

  cout << "expected 3ring events electron/muon=" 
    << expected_3ring_events_electron << "/" << expected_3ring_events_muon << endl;

  for(int r=0;r<r_max;r++){
    for(int mu=0;mu<mu_max;mu++){
      for(int m=0;m<m_max;m++){
        cout << "cut nring/mulike/michel=" << r << "/" << mu << "/" << m << endl;
        cout << "total_with_ntag=" << total_with_ntag[r][mu][m] << endl;
        cout << "lowerSR_without_ntag=" << lowerSR_without_ntag[r][mu][m] << endl;
        cout << "higherSR_without_ntag=" << higherSR_without_ntag[r][mu][m] << endl;
      }
    }
  }


  if(kMakeNtuple){
    if(kDebugMode) cout << "Write ntuple" << endl;
    //TFile *file = new TFile(output_ntuple.c_str(),"recreate");
    otree->Write();
    ofile->Close();
  }

  else m_hSvc.WriteOutput();

}


void OscNtupleManager::Process_pepi(){

  //apply selection here

  if(process_input!="fcmc" && wallv(0)<200) return;

  pass_cut[0][0]=true;
  if ( evis(0) > 30.0 && wall(0) > 200.0 && nhitac(0) <  nhitac_cut[skgen] ){//FC & FV
    pass_cut[1][0]=true;
  }

  //# of cherenkov ring
  if(nring(0)==2) pass_cut[2][0]=true;
  if(nring(0)==3) pass_cut[2][1]=true;
  if(nring(0)==2 || nring(0)==3) pass_cut[2][2]=true;

  //# of mu-like ring
  if(n_mulike_pattern==0) pass_cut[3][0]=true;

  // michel (decay) electron cut
  if(nDecayE==0) pass_cut[4][0]=true;

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

  if(kAllHist) MakeCutFlow();
  else MakeValidationPlot();

}

void OscNtupleManager::Process_pmupi(){

  //apply selection here

  if(process_input!="fcmc" && wallv(0)<200) return;

  pass_cut[0][0]=true;
  if ( evis(0) > 30.0 && wall(0) > 200.0 && nhitac(0) <  nhitac_cut[skgen] ){//FC & FV
    pass_cut[1][0]=true;
  }

  //# of cherenkov ring
  if(nring(0)==2) pass_cut[2][0]=true;
  if(nring(0)==3) pass_cut[2][1]=true;
  if(nring(0)==2 || nring(0)==3) pass_cut[2][2]=true;

  //# of mu-like ring
  if(n_mulike_pattern==1) pass_cut[3][0]=true;

  // michel (decay) electron cut
  if(nDecayE==1) pass_cut[4][0]=true;

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

  if(kAllHist) MakeCutFlow();
  else MakeValidationPlot();

}

void OscNtupleManager::Process_peee(){

  //apply selection here

  if(process_input!="fcmc" && wallv(0)<200) return;

  pass_cut[0][0]=true;
  if ( evis(0) > 30.0 && wall(0) > 200.0 && nhitac(0) <  nhitac_cut[skgen] ){//FC & FV
    pass_cut[1][0]=true;
  }

  //# of cherenkov ring
  if(nring(0)==2) pass_cut[2][0]=true;
  if(nring(0)==3) pass_cut[2][1]=true;
  if(nring(0)==2 || nring(0)==3) pass_cut[2][2]=true;

  //# of mu-like ring
  if(n_mulike_angle==0) pass_cut[3][0]=true;

  // michel (decay) electron cut
  if(nDecayE==0) pass_cut[4][0]=true;

  vector<TLorentzVector> e_cand;
  e_cand.clear();
  for(int r=0;r<nRing;r++)
    if( ( sqrt(fabs(probms(r,1)))-sqrt(fabs(probms(r,2))) ) < 0) e_cand.push_back(GetTLorentzVectorRing(r,1));

  float min_diff=99999999;
  if(e_cand.size()==3){
    for(unsigned int r1=0;r1<e_cand.size()-1;r1++){
      for(unsigned int r2=r1+1;r2<e_cand.size();r2++){
        float mass_pi0_reco = (e_cand[r1] + e_cand[r2]).M();
        float diff = fabs(mass_pi0_reco-135);
        if(diff < min_diff){
          min_diff = diff;
          closest_mass_pi0_reco = mass_pi0_reco;
        }
      }
    }
  }
  //if(nRing==3 && (closest_mass_pi0_reco<85 || closest_mass_pi0_reco>185) ) pass_cut[5][0]=true;
  //if(nRing==2) pass_cut[5][0]=true;

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

  if(process_input!="fcmc" && wallv(0)<200) return;

  pass_cut[0][0]=true;
  if ( evis(0) > 30.0 && wall(0) > 200.0 && nhitac(0) <  nhitac_cut[skgen] ){//FC & FV
    pass_cut[1][0]=true;
  }

  if(nring(0)==3) pass_cut[2][0]=true;

  if(n_mulike_angle==3) pass_cut[3][0]=true;

  // michel (decay) electron cut
  if(nDecayE==1) pass_cut[4][0]=true;
  if(nDecayE==2) pass_cut[4][1]=true;
  if(nDecayE==3) pass_cut[4][2]=true;
  if(nDecayE==2 || nDecayE==3) pass_cut[4][3]=true;
  if(nDecayE==1 || nDecayE==2 || nDecayE==3) pass_cut[4][4]=true;

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

  if (kMakeNtuple) MakeNtuple();
  if(kAllHist) MakeCutFlow();
  else MakeValidationPlot();

}

void OscNtupleManager::Process_pemumu(){

  //apply selection here

  if(process_input!="fcmc" && wallv(0)<200) return;

  pass_cut[0][0]=true;
  if ( evis(0) > 30.0 && wall(0) > 200.0 && nhitac(0) <  nhitac_cut[skgen] ){//FC & FV
    pass_cut[1][0]=true;
  }

  //# of cherenkov ring
  if(nring(0)==2) pass_cut[2][0]=true;
  if(nring(0)==3) pass_cut[2][1]=true;
  if(nring(0)==2 || nring(0)==3) pass_cut[2][2]=true;

  //#of mu-like ring
  if(n_mulike_angle==2) pass_cut[3][0]=true;

  // michel (decay) electron cut
  if(nDecayE==1) pass_cut[4][0]=true;
  if(nDecayE==2) pass_cut[4][1]=true;
  if(nDecayE==1 || nDecayE==2) pass_cut[4][2]=true;

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

  if(process_input!="fcmc" && wallv(0)<200) return;

  pass_cut[0][0]=true;
  if ( evis(0) > 30.0 && wall(0) > 200.0 && nhitac(0) <  nhitac_cut[skgen] ){//FC & FV
    pass_cut[1][0]=true;
  }

  //# of cherenkov ring
  if(nring(0)==2) pass_cut[2][0]=true;
  if(nring(0)==3) pass_cut[2][1]=true;
  if(nring(0)==2 || nring(0)==3) pass_cut[2][2]=true;

  //#of mu-like ring
  if(n_mulike_angle==1) pass_cut[3][0]=true;

  // michel (decay) electron cut
  if(nDecayE==0) pass_cut[4][0]=true;
  if(nDecayE==1) pass_cut[4][1]=true;
  if(nDecayE==0 || nDecayE==1) pass_cut[4][2]=true;

  vector<TLorentzVector> e_cand;
  TLorentzVector mu_cand;
  e_cand.clear();
  for(int r=0;r<nRing;r++){
    if( ( sqrt(fabs(probms(r,1)))-sqrt(fabs(probms(r,2))) ) < 0) e_cand.push_back(GetTLorentzVectorRing(r,1));
    else mu_cand = GetTLorentzVectorRing(r,2);
  }

  TLorentzVector total_vec,two_elike_vec;
  if(n_elike_angle==2) {
    total_vec = e_cand[0] + e_cand[1] + mu_cand;
    two_elike_vec = e_cand[0] + e_cand[1];
  }
  if(n_elike_angle==1) total_vec = e_cand[0] + mu_cand;
  total_mass = total_vec.M();
  total_mom = total_vec.P();
  closest_mass_pi0_reco = two_elike_vec.M();
  
  if(nRing==3 && (closest_mass_pi0_reco<85 || closest_mass_pi0_reco>185) ) pass_cut[5][0]=true;
  if(nRing==2) pass_cut[5][0]=true;

  //pass_cut[5][0]=true;//no pi0 mass cut


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

void OscNtupleManager::Process_single(){

  m_hSvc.h1D("ring_counting_likelihood","","")->Fill(Dlfct(0));
  m_hSvc.h1D("true_mom_lepton","","")->Fill(pmomv(0));
  m_hSvc.h1D("true_dirx_lepton","","")->Fill(dirv(0,0));
  m_hSvc.h1D("true_diry_lepton","","")->Fill(pmomv(0,1));
  m_hSvc.h1D("true_dirz_lepton","","")->Fill(pmomv(0,2));
  m_hSvc.h1D("nRing","","")->Fill(nRing);
  m_hSvc.h1D("nMulikeRing_angle","","")->Fill(n_mulike_angle);

  TLorentzVector this_lepton = GetTLorentzVectorVector(0);
  float closest_angle=9999,closest_ring_prob_angle=9999;
  int closest_ring_id=-1;
  for(int r=0;r<nRing;r++){
    float prob_angle = sqrt(fabs(probms(r,1)))-sqrt(fabs(probms(r,2)));
    int pid_angle = (prob_angle<0)? 0 : 2;//e-like or mu-like
    TLorentzVector this_ring = GetTLorentzVectorRing(r,pid_angle);
    float angle_lep_ring = this_ring.Angle(this_lepton.Vect())*180./3.14159;
    m_hSvc.h1D("true_angle_lepton_and_ring","","")->Fill(angle_lep_ring);
    if(angle_lep_ring<20) {
      if(angle_lep_ring<closest_angle){
        closest_angle = angle_lep_ring;
        closest_ring_id = r;
        closest_ring_prob_angle = prob_angle;
      }
    }
  }

  if(kDebugMode) cout << "closest_ring_id=" << closest_ring_id << endl;
  if(ipv(0)==2 || ipv(0)==3) m_hSvc.h1D("true_mom_electron","","")->Fill(pmomv(0));
  if(ipv(0)==5 || ipv(0)==6) m_hSvc.h1D("true_mom_muon","","")->Fill(pmomv(0));
  if(closest_ring_id!=-1){//matching success !!
    float prob_angle = sqrt(fabs(probms(closest_ring_id,1)))-sqrt(fabs(probms(closest_ring_id,2)));
    float residual_emom = (amome(closest_ring_id)-pmomv(0))/pmomv(0);
    float residual_mmom = (amomm(closest_ring_id)-pmomv(0))/pmomv(0);
    if(ipv(0)==2 || ipv(0)==3){//electron
      m_hSvc.h1D("true_mom_electron_match_ring","","")->Fill(pmomv(0));
      if(prob_angle<0) m_hSvc.h1D("true_mom_electron_match_ring_angle_elike","","")->Fill(pmomv(0));
    }
    if(ipv(0)==5 || ipv(0)==6){//muon
      m_hSvc.h1D("true_mom_muon_match_ring","","")->Fill(pmomv(0));
      if(prob_angle>0) m_hSvc.h1D("true_mom_muon_match_ring_angle_mulike","","")->Fill(pmomv(0));
    }
    for(int p=0;p<6;p++){
      if(pmomv(0)>100*p && pmomv(0)<100+100*p){
        if(ipv(0)==2 || ipv(0)==3) {//electron
          m_hSvc.h1D(Form("residual_emom_mom%d_%d",100*p,100+100*p),"","")->Fill(residual_emom);
          m_hSvc.h1D(Form("prob_angle_electron_mom%d_%d",100*p,100+100*p),"","")->Fill(prob_angle);
        }
        if(ipv(0)==5 || ipv(0)==6) {//muon
          m_hSvc.h1D(Form("residual_mmom_mom%d_%d",100*p,100+100*p),"","")->Fill(residual_mmom);
          m_hSvc.h1D(Form("prob_angle_muon_mom%d_%d",100*p,100+100*p),"","")->Fill(prob_angle);
        }
      }
    }
  }

  if (kMakeNtuple) MakeNtuple();

}

void OscNtupleManager::Process_subgev_multiring(){

  if(kDebugMode) cout << "Process_subgev_multiring" << endl;
  
  if ( evis(0) < 30.0 || evis(0) > 1330.0 || nhitac(0) >  nhitac_cut[skgen] || nring(0)<2) 
    return;//FC subgev multi ring w/o dwall cut

  int dwall_thr[3] = {50,100,150};
  for(int d=0;d<3;d++){
    if(wall(0)>dwall_thr[d]){
      m_hSvc.h1D(Form("distance_to_wall_thr%d",dwall_thr[d]),"","")->Fill(wall(0),weight*osc_weight);
      m_hSvc.h1D(Form("distance_to_wall_thr%d_nring%d",dwall_thr[d],nRing),"","")->Fill(wall(0),weight*osc_weight);
    }
  }

  //w/ Dwall cut for decayE systematic in multi ring sample
  if(wall(0)>200){
    m_hSvc.h1D(Form("n_michel_electron_mulike%d",n_mulike_angle),"","")->Fill(nDecayE,weight*osc_weight);
    int type=-1;
    if(abs(ipnu(0))==12) type=0;//nue or nuebar
    if(abs(ipnu(0))==14) type=1;//numu or numubar
    if(mode(0)>0) m_hSvc.h1D(Form("n_michel_electron_mulike%d_type%d_mode_pos%d",n_mulike_angle,type,abs(mode(0))),"","")->Fill(nDecayE,weight*osc_weight);
    else m_hSvc.h1D(Form("n_michel_electron_mulike%d_type%d_mode_neg%d",n_mulike_angle,type,abs(mode(0))),"","")->Fill(nDecayE,weight*osc_weight);
  }

}

void OscNtupleManager::Process_subgev_onemulike(){

  if(kDebugMode) cout << "Process_subgev_onemulike" << endl;

  if ( evis(0) < 30.0 || evis(0) > 1330.0 || nhitac(0) > nhitac_cut[skgen] 
      || wall(0) < 200 || nring(0) != 1 || n_mulike_angle != 1) 
    return;//FC subgev one mu-like ring 

    m_hSvc.h1D("n_michel_electron","","")->Fill(nDecayE,weight*osc_weight);
    m_hSvc.h1D("nRing","","")->Fill(nRing,weight*osc_weight);
    m_hSvc.h1D("nMulikeRing_angle","","")->Fill(n_mulike_angle,weight*osc_weight);
    int type=-1;
    if(abs(ipnu(0))==12) type=0;//nue or nuebar
    if(abs(ipnu(0))==14) type=1;//numu or numubar
    if(mode(0)>0) m_hSvc.h1D(Form("n_michel_electron_type%d_mode_pos%d",type,abs(mode(0))),"","")->Fill(nDecayE,weight*osc_weight);
    else m_hSvc.h1D(Form("n_michel_electron_type%d_mode_neg%d",type,abs(mode(0))),"","")->Fill(nDecayE,weight*osc_weight);

}

void OscNtupleManager::Process_subgev_oneelike(){

  if(kDebugMode) cout << "Process_subgev_oneelike" << endl;

  if ( evis(0) < 30.0 || evis(0) > 1330.0 || nhitac(0) > nhitac_cut[skgen] 
      || wall(0) < 200 || nring(0) != 1 || n_elike_angle != 1) 
    return;//FC subgev one mu-like ring 

    m_hSvc.h1D("n_michel_electron","","")->Fill(nDecayE,weight*osc_weight);
    m_hSvc.h1D("nRing","","")->Fill(nRing,weight*osc_weight);
    m_hSvc.h1D("nMulikeRing_angle","","")->Fill(n_mulike_angle,weight*osc_weight);

}

void OscNtupleManager::ZeroStructure()
{
  if(kDebugMode) cout << "OscNtupleManager::ZeroStructure" << endl;

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
  if(process_input=="single_mu_def" || process_input=="single_mu_miura") nNeutron = 0;
  else if(skgen == SK4) nNeutron = ntag_nn(0);//only for sk4
  else nNeutron = 0;

  //find the events with captured neutron
  n_true_neutron=0;
  if(process_input=="fcmc"){
    //cout << "# of secondary particle is " << nscndprt(0) << endl;
    for(int n=0;n<nscndprt(0);n++){
      //cout << "#/pid/ind=" << n << "/" << iprtscnd(n) << "/" << lmecscnd(n) << endl;
      if(iprtscnd(n)==100045 && lmecscnd(n)==18) n_true_neutron++;
    }
  }
  for(int n=0;n<11;n++) n_tagged_neutron_exp[n] = 0;
  TRandom *generator = new TRandom();
  generator->SetSeed(0);
  for(int e=0;e<11;e++){
    float efficiency = 0.1*e;
    for(int n=0;n<n_true_neutron;n++){
      float rnd = generator->Rndm();
      if(rnd<efficiency) n_tagged_neutron_exp[e]++;
    }
  }

  for(int c=0;c<10;c++) 
    for(int c2=0;c2<5;c2++)
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
    cout << " Error - unsupported mode : " << mode  << " exiting " << endl;
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
    //cout << " MER ID: " << mer_id << " " << llelectron << endl;

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
  //cout << "ipnu/NuOscillatedTo=" << ipnu(0) << "/" << NuOscillatedTo << endl;
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
  //cout << "this_itype=" << this_itype << endl;
  //cout << "skgen/EndOfTypes=" << skgen << "/" << EndOfTypes << endl;

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
  //cout << "FactorE/FactorMu=" << FactorE << "/" << FactorMu << endl;

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

  //cout << "FullPathAve=" << FullPathAve << endl;

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
  //cout << "Type=" << Type << endl;

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
  //cout << "OscNtupleManager::GetMCweight" << endl;

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
  //cout << "skgen/weightx=" << skgen << "/" << weightx << endl;

  return weightx;

}

double OscNtupleManager::GetHondaFluxRatio( int NuType )
{
  //cout << "GetHondaFluxRatio" << endl;
  float Solar    = 0.5;
  float d_pnu = pnu(0);
  float d_dirnu[3];
  for( int i = 0 ; i < 3 ; i++ ){
    d_dirnu[i] = dirnu(0,i);
    //cout << "d_dirnu=" << d_dirnu[i] << endl;
  }

  int   type[] = {12, 14};
  double flxho[2];
  for( int i = 0 ; i < 2 ; i++ )
  {
    int nu = type[i] * ( ipnu(0) < 1 ? -1 : 1 );
    //cout << "d_pnu/d_dir_nu/Solar/nu=" << d_pnu << "/" << d_dirnu << "/" << Solar << "/" << nu << endl;
    //flxho[i] = fnhonfx11_( &d_pnu, d_dirnu, &Solar, &nu ) ;
    flxho[i] = m_fort->calc_flux( d_pnu, d_dirnu, Solar, nu ) ;
    //cout << "flxho=" << flxho[i] << endl;
  }

  if ( NuType == 1 )   
    return (double) flxho[1] / flxho[0];

  if ( NuType == 2 || NuType == 3 )
    return (double) flxho[1] / flxho[0];

  std::cerr << "Returning 0 from SKEventParser::GetHondaFluxRatio for type " << NuType << endl;
  return 0;

}

void OscNtupleManager::BuildEnergyFriend(int type){
  //cout << "OscNtupleManager::BuildEnergyFriend" << endl;
  int zbin = GetZBin(-1.*dir(2,0),type);
  //cout << "zbin=" << zbin << endl;


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
  if(kDebugMode) cout << "FillNtuple" << endl;

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
    cout << "interacion_type=" << interaction_type << endl;
    cout << "weight mc/osc=" << mc_weight << "/" << osc_weight << endl;
  }
  if(interaction_type!=-1)  {
    m_hSvc.h1D(Form("nRing_type%d",interaction_type),"","")->Fill(nring(0),weight);
    m_hSvc.h1D(Form("nRing_type%d_osc",interaction_type),"","")->Fill(nring(0),weight*osc_weight);
  }
  if(interaction_type>=1 && interaction_type<=10 && interaction_type!=7){//single ring FC event
    //cout << "single ring event?" << endl;
    //cout << "nring=" << nring(0) << endl;
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
  
  if(kOutsideSR){//reject events inside the signal region
    if(total_mass>800 && total_mass<1050 && total_mom<250) return;
  }

  for(int r=0;r<r_max;r++){// # of ring
    for(int mu=0;mu<mu_max;mu++){// # of mu-like ring
      for(int m=0;m<m_max;m++){//# of michel electron
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

void OscNtupleManager::MakeBasicPlot(int c, int r, int mu, int p){//cut #, michel e cut, 

  if(process_input!="fcdt" || total_mass<800 || total_mass>1050 || total_mom>250){//for blind analysis
    m_hSvc.h1D(Form("cut_flow_nring%d_mulike%d_michel%d",r,mu,p),"","")->Fill(c,weight*osc_weight);
    m_hSvc.h1D(Form("cut_flow_nring%d_mulike%d_michel%d_fp%d",r,mu,p,is_free_proton),"","")->Fill(c,weight*osc_weight);
  }
  for(int e=0;e<11;e++){
    int eff = 10*e;
    m_hSvc.h1D(Form("n_tagged_neutron_exp_eff%d_cut%d_nring%d_mulike%d_michel%d",eff,c,r,mu,p),"","")->Fill(n_tagged_neutron_exp[e],weight*osc_weight);
  }
  m_hSvc.h1D(Form("interaction_mode_cut%d_nring%d_mulike%d_michel%d",c,r,mu,p),"","")->Fill(true_mode,weight*osc_weight);
  m_hSvc.h1D(Form("n_true_neutron_cut%d_nring%d_mulike%d_michel%d",c,r,mu,p),"","")->Fill(n_true_neutron,weight*osc_weight);
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
  m_hSvc.h1D(Form("mass_pi0_reco_cut%d_nring%d_mulike%d_michel%d",c,r,mu,p),"","")->Fill(closest_mass_pi0_reco,weight*osc_weight);

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
  m_hSvc.h1D(Form("mass_pi0_reco_cut%d_nring%d_mulike%d_michel%d_fp%d",c,r,mu,p,is_free_proton),"","")->Fill(closest_mass_pi0_reco,weight*osc_weight);
  if(process_input!="fcdt" || total_mass<800 || total_mass>1050){//for blind analysis
    m_hSvc.h1D(Form("mass_proton_reco_cut%d_nring%d_mulike%d_michel%d",c,r,mu,p),"","")->Fill(total_mass,weight*osc_weight);
    m_hSvc.h1D(Form("mass_proton_reco_cut%d_nring%d_mulike%d_michel%d_fp%d",c,r,mu,p,is_free_proton),"","")->Fill(total_mass,weight*osc_weight);
    m_hSvc.h1D(Form("mass_two_elike_reco_cut%d_nring%d_mulike%d_michel%d",c,r,mu,p),"","")->Fill(closest_mass_pi0_reco,weight*osc_weight);
    m_hSvc.h1D(Form("mass_two_elike_reco_cut%d_nring%d_mulike%d_michel%d_fp%d",c,r,mu,p,is_free_proton),"","")->Fill(closest_mass_pi0_reco,weight*osc_weight);
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
    //cout << "Fill!!" << endl;
    int type=-1;
    if(abs(ipnu(0))==12) type=0;//nue or nuebar
    if(abs(ipnu(0))==14) type=1;//numu or numubar
    if(mode(0)>0){
      m_hSvc.h1D(Form("cut_flow_nring%d_mulike%d_michel%d_type%d_mode_pos%d",r,mu,p,type,abs(mode(0))),"","")->Fill(c,weight*osc_weight);
      m_hSvc.h1D(Form("nRing_cut%d_nring%d_mulike%d_michel%d_type%d_mode_pos%d",c,r,mu,p,type,abs(mode(0))),"","")->Fill(nRing,weight*osc_weight);
      m_hSvc.h1D(Form("nElikeRing_angle_nring%d_cut%d_nring%d_mulike%d_michel%d_type%d_mode_pos%d",nRing,c,r,mu,p,type,abs(mode(0))),"","")->Fill(n_elike_angle,weight*osc_weight);
      m_hSvc.h1D(Form("nMulikeRing_angle_nring%d_cut%d_nring%d_mulike%d_michel%d_type%d_mode_pos%d",nRing,c,r,mu,p,type,abs(mode(0))),"","")->Fill(n_mulike_angle,weight*osc_weight);
      m_hSvc.h1D(Form("n_michel_electron_cut%d_nring%d_mulike%d_michel%d_type%d_mode_pos%d",c,r,mu,p,type,abs(mode(0))),"","")->Fill(nDecayE,weight*osc_weight);
      m_hSvc.h1D(Form("ntag_multiplicity_cut%d_nring%d_mulike%d_michel%d_type%d_mode_pos%d",c,r,mu,p,type,abs(mode(0))),"","")->Fill(nNeutron,weight*osc_weight);
      m_hSvc.h1D(Form("mass_proton_reco_cut%d_nring%d_mulike%d_michel%d_type%d_mode_pos%d",c,r,mu,p,type,abs(mode(0))),"","")->Fill(total_mass,weight*osc_weight);
      m_hSvc.h1D(Form("mass_two_elike_reco_cut%d_nring%d_mulike%d_michel%d_type%d_mode_pos%d",c,r,mu,p,type,abs(mode(0))),"","")->Fill(closest_mass_pi0_reco,weight*osc_weight);
      m_hSvc.h1D(Form("mom_proton_reco_cut%d_nring%d_mulike%d_michel%d_type%d_mode_pos%d",c,r,mu,p,type,abs(mode(0))),"","")->Fill(total_mom,weight*osc_weight);
    }
    else {//negative mode number
      m_hSvc.h1D(Form("cut_flow_nring%d_mulike%d_michel%d_type%d_mode_neg%d",r,mu,p,type,abs(mode(0))),"","")->Fill(c,weight*osc_weight);
      m_hSvc.h1D(Form("nRing_cut%d_nring%d_mulike%d_michel%d_type%d_mode_neg%d",c,r,mu,p,type,abs(mode(0))),"","")->Fill(nRing,weight*osc_weight);
      m_hSvc.h1D(Form("nRing_cut%d_nring%d_mulike%d_michel%d_type%d_mode_neg%d",c,r,mu,p,type,abs(mode(0))),"","")->Fill(nRing,weight*osc_weight);
      m_hSvc.h1D(Form("nElikeRing_angle_nring%d_cut%d_nring%d_mulike%d_michel%d_type%d_mode_neg%d",nRing,c,r,mu,p,type,abs(mode(0))),"","")->Fill(n_elike_angle,weight*osc_weight);
      m_hSvc.h1D(Form("nMulikeRing_angle_nring%d_cut%d_nring%d_mulike%d_michel%d_type%d_mode_neg%d",nRing,c,r,mu,p,type,abs(mode(0))),"","")->Fill(n_mulike_angle,weight*osc_weight);
      m_hSvc.h1D(Form("n_michel_electron_cut%d_nring%d_mulike%d_michel%d_type%d_mode_neg%d",c,r,mu,p,type,abs(mode(0))),"","")->Fill(nDecayE,weight*osc_weight);
      m_hSvc.h1D(Form("ntag_multiplicity_cut%d_nring%d_mulike%d_michel%d_type%d_mode_neg%d",c,r,mu,p,type,abs(mode(0))),"","")->Fill(nNeutron,weight*osc_weight);
      m_hSvc.h1D(Form("mass_proton_reco_cut%d_nring%d_mulike%d_michel%d_type%d_mode_neg%d",c,r,mu,p,type,abs(mode(0))),"","")->Fill(total_mass,weight*osc_weight);
      m_hSvc.h1D(Form("mass_two_elike_reco_cut%d_nring%d_mulike%d_michel%d_type%d_mode_neg%d",c,r,mu,p,type,abs(mode(0))),"","")->Fill(closest_mass_pi0_reco,weight*osc_weight);
      m_hSvc.h1D(Form("mom_proton_reco_cut%d_nring%d_mulike%d_michel%d_type%d_mode_neg%d",c,r,mu,p,type,abs(mode(0))),"","")->Fill(total_mom,weight*osc_weight);
    }
  }
  graph_point_2++;

}

void OscNtupleManager::MakeValidationPlot(){
  if(kDebugMode) cout << "OscNtupleManager::MakeValidationPlot" << endl;

  float initial_pmom = 0.,initial_pdirx=0.,initial_pdiry=0.,initial_pdirz=0.;
  if(process_input=="fcmc"){
    if(ipv(1)==14) {
      initial_pmom = pmomv(1);
      initial_pdirx = dirv(1,0);
      initial_pdiry = dirv(1,1);
      initial_pdirz = dirv(1,2);
    }
  }
  else {//for proton decay events
    initial_pmom = pmomv(0);
    initial_pdirx = dirv(0,0);
    initial_pdiry = dirv(0,1);
    initial_pdirz = dirv(0,2);
  }
  if(process_input=="fcmc"){
    if(initial_pmom>1e-2) {
      m_hSvc.h1D("fermi_momentum","","")->Fill(initial_pmom,weight);
      m_hSvc.h1D("fermi_dirx","","")->Fill(initial_pdirx,weight);
      m_hSvc.h1D("fermi_diry","","")->Fill(initial_pdiry,weight);
      m_hSvc.h1D("fermi_dirz","","")->Fill(initial_pdirz,weight);
    }
  } else if (hstate(0)!=0 && hstate(0)!=4){
    m_hSvc.h1D("fermi_momentum","","")->Fill(initial_pmom,weight);
    m_hSvc.h1D("fermi_dirx","","")->Fill(initial_pdirx,weight);
    m_hSvc.h1D("fermi_diry","","")->Fill(initial_pdiry,weight);
    m_hSvc.h1D("fermi_dirz","","")->Fill(initial_pdirz,weight);
  }


  //Basic plot
  if(process_input!="fcmc" && process_input!="fcdt") m_hSvc.h1D("hstate","","")->Fill(hstate(0),weight);

  //truth&ring matching
  //ordering by momentum 
  float min_mom=99999999, mid_mom=99999999, max_mom=99999999;
  int min_mom_id=-1,mid_mom_id=-1,max_mom_id=-1;
  int matched_ring_id[4];
  float matched_ring_prob_angle[4];
  TLorentzVector total_true_vec;
  for(int v=1;v<4;v++){//loop for 3 leptons
    if(kDebugMode) cout << "lepton mom/pid/x/y/z=" << pmomv(v) << "/" << ipv(v) << "/" << dirv(v,0) << "/" << dirv(v,1) << "/" << dirv(v,2) << endl;
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
    float closest_angle=9999,closest_ring_prob_angle=9999;
    int closest_ring_id=-1;
    for(int r=0;r<nRing;r++){
      float prob_angle = sqrt(fabs(probms(r,1)))-sqrt(fabs(probms(r,2)));
      int pid_angle = (prob_angle<0)? 0 : 2;//e-like or mu-like
      TLorentzVector this_ring = GetTLorentzVectorRing(r,pid_angle);
      float angle_lep_ring = this_ring.Angle(this_lepton.Vect())*180./3.14159;
      m_hSvc.h1D(Form("true_angle_lepton_and_ring_nring%d",nRing),"","")->Fill(angle_lep_ring);
      if(angle_lep_ring<20) {
        if(angle_lep_ring<closest_angle){
          closest_angle = angle_lep_ring;
          closest_ring_id = r;
          closest_ring_prob_angle = prob_angle;
        }
      }
    }
    matched_ring_id[v] = closest_ring_id;
    matched_ring_prob_angle[v] = closest_ring_prob_angle;
    if(kDebugMode) cout << "closest_ring_id=" << closest_ring_id << endl;
  }
  if(kDebugMode) cout << "minmom id/mom=" << min_mom_id << "/" << min_mom << endl;
  if(kDebugMode) cout << "midmom id/mom=" << mid_mom_id << "/" << mid_mom << endl;
  if(kDebugMode) cout << "maxmom id/mom=" << max_mom_id << "/" << max_mom << endl;
  TLorentzVector min_mom_lepton = GetTLorentzVectorVector(min_mom_id);
  TLorentzVector mid_mom_lepton = GetTLorentzVectorVector(mid_mom_id);
  TLorentzVector max_mom_lepton = GetTLorentzVectorVector(max_mom_id);
  float true_angle_min_mid_lepton = min_mom_lepton.Angle(mid_mom_lepton.Vect())*180./3.14159;
  float true_angle_min_max_lepton = min_mom_lepton.Angle(max_mom_lepton.Vect())*180./3.14159;
  float true_angle_mid_max_lepton = mid_mom_lepton.Angle(max_mom_lepton.Vect())*180./3.14159;
  m_hSvc.h1D("true_angle_min_mid_lepton","","")->Fill(true_angle_min_mid_lepton,weight);
  m_hSvc.h1D("true_angle_min_max_lepton","","")->Fill(true_angle_min_max_lepton,weight);
  m_hSvc.h1D("true_angle_mid_max_lepton","","")->Fill(true_angle_mid_max_lepton,weight);

  //make histgram
  m_hSvc.h1D("ring_counting_likelihood","","")->Fill(Dlfct(0));
  m_hSvc.h1D(Form("true_min_mom_lepton_fp%d",is_free_proton),"","")->Fill(min_mom,weight);
  m_hSvc.h1D(Form("true_mid_mom_lepton_fp%d",is_free_proton),"","")->Fill(mid_mom,weight);
  m_hSvc.h1D(Form("true_max_mom_lepton_fp%d",is_free_proton),"","")->Fill(max_mom,weight);
  for(int c=1;c<4;c++){
    m_hSvc.h1D(Form("true_mom_lepton_fp%d",is_free_proton),"","")->Fill(pmomv(c),weight);
    m_hSvc.h1D(Form("true_mom_lepton_nring%d",nRing),"","")->Fill(pmomv(c));
    if(ipv(c)==2 || ipv(c)==3){
      m_hSvc.h1D("true_mom_electron","","")->Fill(pmomv(c));
      m_hSvc.h1D("prob_angle_electron","","")->Fill(matched_ring_prob_angle[c]);
    }
    if(ipv(c)==5 || ipv(c)==6) {
      m_hSvc.h1D("true_mom_muon","","")->Fill(pmomv(c));
      m_hSvc.h1D(Form("true_mom_muon_nring%d",nRing),"","")->Fill(pmomv(c));
      m_hSvc.h1D(Form("true_mom_muon_nring%d_mulike%d",nRing,n_mulike_angle),"","")->Fill(pmomv(c));
      m_hSvc.h1D("prob_angle_muon","","")->Fill(matched_ring_prob_angle[c]);
    }
    if(matched_ring_id[c]){//matching success !!
      float residual_emom = (amome(matched_ring_id[c])-pmomv(c))/pmomv(c);
      float residual_mmom = (amomm(matched_ring_id[c])-pmomv(c))/pmomv(c);
      if(ipv(c)==2 || ipv(c)==3){//electron
        m_hSvc.h1D(Form("true_mom_electron_match_ring_nring%d",nRing),"","")->Fill(pmomv(c));
        if(matched_ring_prob_angle[c]<0) m_hSvc.h1D(Form("true_mom_electron_match_ring_angle_elike_nring%d",nRing),"","")->Fill(pmomv(c));
      }
      if(ipv(c)==5 || ipv(c)==6){//muon
        m_hSvc.h1D(Form("true_mom_muon_match_ring_nring%d",nRing),"","")->Fill(pmomv(c));
        if(matched_ring_prob_angle[c]>0) m_hSvc.h1D(Form("true_mom_muon_match_ring_angle_mulike_nring%d",nRing),"","")->Fill(pmomv(c));
      }
      for(int p=0;p<6;p++){
        if(pmomv(c)>100*p && pmomv(c)<100+100*p){
          if(ipv(c)==2 || ipv(c)==3) {//electron
            m_hSvc.h1D(Form("residual_emom_mom%d_%d_nring%d",100*p,100+100*p,nRing),"","")->Fill(residual_emom);
            m_hSvc.h1D(Form("prob_angle_electron_mom%d_%d",100*p,100+100*p),"","")->Fill(matched_ring_prob_angle[c]);
          }
          if(ipv(c)==5 || ipv(c)==6) {//muon
            m_hSvc.h1D(Form("residual_mmom_mom%d_%d_nring%d",100*p,100+100*p,nRing),"","")->Fill(residual_mmom);
            m_hSvc.h1D(Form("prob_angle_muon_mom%d_%d",100*p,100+100*p),"","")->Fill(matched_ring_prob_angle[c]);
          }
        }
      }
    }
  }

  float total_true_mass = total_true_vec.M();
  float total_true_mom = total_true_vec.P();
  float angle_min_mom_lepton_proton = min_mom_lepton.Angle(total_true_vec.Vect())*180./3.14159;
  float angle_mid_mom_lepton_proton = mid_mom_lepton.Angle(total_true_vec.Vect())*180./3.14159;
  float angle_max_mom_lepton_proton = max_mom_lepton.Angle(total_true_vec.Vect())*180./3.14159;
  float residual_total_mass = (total_mass-total_true_mass)/total_true_mass;
  float diff_total_mass = total_mass-total_true_mass;
  float diff_total_mom = total_mom-total_true_mom;
  m_hSvc.h1D(Form("angle_min_mom_lepton_proton_fp%d",is_free_proton),"","")->Fill(angle_min_mom_lepton_proton,weight);
  m_hSvc.h1D(Form("angle_mid_mom_lepton_proton_fp%d",is_free_proton),"","")->Fill(angle_mid_mom_lepton_proton,weight);
  m_hSvc.h1D(Form("angle_max_mom_lepton_proton_fp%d",is_free_proton),"","")->Fill(angle_max_mom_lepton_proton,weight);
  m_hSvc.h1D(Form("total_true_mass_fp%d",is_free_proton),"","")->Fill(total_true_mass,weight);
  m_hSvc.h1D(Form("total_true_mom_fp%d",is_free_proton),"","")->Fill(total_true_mom,weight);
  m_hSvc.h1D(Form("total_mass_nring%d_mulike%d_fp%d",nRing,n_mulike_angle,is_free_proton),"","")->Fill(total_mass,weight);
  m_hSvc.h1D(Form("total_mom_nring%d_mulike%d_fp%d",nRing,n_mulike_angle,is_free_proton),"","")->Fill(total_mom,weight);
  m_hSvc.h1D(Form("residual_total_mass_nring%d_mulike%d_fp%d",nRing,n_mulike_angle,is_free_proton),"","")->Fill(residual_total_mass);
  m_hSvc.h1D(Form("diff_total_mass_nring%d_mulike%d_fp%d",nRing,n_mulike_angle,is_free_proton),"","")->Fill(diff_total_mass);
  m_hSvc.h1D(Form("diff_total_mom_nring%d_mulike%d_fp%d",nRing,n_mulike_angle,is_free_proton),"","")->Fill(diff_total_mom);

  float vertex_r_ring = sqrt(pos(0)*pos(0)+pos(1)*pos(1)+pos(2)*pos(2));
  float vertex_r_true = sqrt(posv(0)*posv(0)+posv(1)*posv(1)+posv(2)*posv(2));
  float diff_vertex_x = pos(0) - posv(0);
  float diff_vertex_y = pos(1) - posv(1);
  float diff_vertex_z = pos(2) - posv(2);
  float diff_vertex_r = vertex_r_ring - vertex_r_true;
  m_hSvc.h1D("diff_vertex_x","","")->Fill(diff_vertex_x);
  m_hSvc.h1D("diff_vertex_y","","")->Fill(diff_vertex_y);
  m_hSvc.h1D("diff_vertex_z","","")->Fill(diff_vertex_z);
  m_hSvc.h1D("diff_vertex_r","","")->Fill(diff_vertex_r);

}

void OscNtupleManager::MakeNtuple(){
  if(kDebugMode) {
    cout << "MakeNtuple" << endl;
    cout << "nRing/nMulike=" << nRing << "/" << n_mulike_angle << endl;
  }
  o_itype = 0 ;  
  o_nring = nRing;
  o_nmulike = n_mulike_angle;
  o_total_mass = total_mass;
  o_total_mom = total_mom;
  o_dlfct = Dlfct(0);
  o_vertex_x = pos(0); 
  o_vertex_y = pos(1); 
  o_vertex_z = pos(2); 
  float min_mom=99999999, mid_mom=99999999, max_mom=99999999;
  for(int r=0;r<nRing;r++){
    //cout << "ring" << r << " amomm=" << amomm(r) << endl;
    float prob_angle = sqrt(fabs(probms(r,1)))-sqrt(fabs(probms(r,2)));
    float prob_pattern = prmslg(r,1) - prmslg(r,2);
    o_prob_angle[r] = prob_angle;
    o_probms_e[r] = probms(r,1);
    o_probms_mu[r] = probms(r,2);
    o_prob_pattern[r] = prob_pattern;
    o_prmslg_e[r] = prmslg(r,1);
    o_prmslg_mu[r] = prmslg(r,2);
    o_mmom[r] = amomm(r);
    o_dir_x[r] = dir(r,0);
    o_dir_y[r] = dir(r,1);
    o_dir_z[r] = dir(r,2);
    o_ang[r] = ang(r);
    o_ange[r] = ange(r);
    o_angm[r] = angm(r);
    if(amomm(r)<min_mom) {
      max_mom=mid_mom;
      mid_mom=min_mom;
      min_mom=amomm(r);
    }
    else if(amomm(r)<mid_mom) {
      max_mom=mid_mom;
      mid_mom=amomm(r);
    }
    else if(amomm(r)<max_mom) {
      max_mom=amomm(r);
    }
  }
  //cout << "amomm min/mid/max=" << min_mom << "/" << mid_mom << "/" << max_mom << endl;
  o_mmom_min = min_mom;
  o_mmom_mid = mid_mom;
  o_mmom_max = max_mom;

  otree->Fill();
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

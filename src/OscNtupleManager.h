#ifndef _OscNtupleManager_
#define _OscNtupleManager_

#include "TLorentzVector.h"
#include "TVector3.h"

#include <string>
#include <vector>

#include "TTree.h"
#include "TFile.h"

#include "tools/DataManager.h"
#include "tools/HistogramService.h"

#include "ESL/pcstopLikelihood.h"
#include "ESL/numu1RLikelihood.h"
#include "ESL/numuLikelihood.h"
#include "ESL/nuebarLikelihood.h"
#include "ESL/mgmreLikelihood.h"
#include "ESL/Pi0Likelihood.h"
#include "src/outputOscStructure.h"
#include "ESL/LikelihoodManager.h"

//#include "osc_types_oct11.h"
#include "src/osc_types.h"

#include "san.sedai/ThreeFlavorOscillator.h"
#include "Prob3++/BargerPropagator.h"

#include "fortran/fort_func.h"

#include "TArrayF.h"

class OscNtupleManager
{

  public:
    // Build Likelihood functions using default,
    // official likelihood*.root files
    OscNtupleManager(DataManager *, int , const char *);

    // Specify your own likelihood*.root files through,
    // card reader
    OscNtupleManager(DataManager *, CardReader *, int , const char *);
   ~OscNtupleManager( );

   HistogramService m_hSvc;
    void LoadWrappers();
    void CreateHist();
    void CreateBranch();
    //std::vector<TH1*> m_TH1_list;

    void SetInputTree( TTree * lt ) { _ltree = lt; }
    void SetOutputBranches();
    void SetOutputFile(TFile* file);
    void SetOutput();
    void SetOutput_ERecon();
    void SetBranch_ERecon();
    void Make_Ntuple();
    void MakeOscillationPlot();
    void MakeBasicPlot(int cut, int nring, int nmulike, int nmichel);
    void MakeCutFlow();
    void MakeCutFlowValidate();//for epi & mupi
    void MakeValidationPlot();//mainly truth plot
    
    void Process(int seed=0, int blossom =-1);
    void SetMode( std::string mode ) { process_mode = mode ; }
    void SetInput( std::string input ) { process_input = input ; }
    void SetOutputNtuple( std::string output ) { output_ntuple = output ; }
    void SetOutputHist( std::string output ) { output_hist = output ; }
    void SetSKX( int sk ) {skgen = sk - 1;}
    void SetLiveTimeWeight(float w) {live_time_weight = w;}

    bool SetEventType( std::string mode , int &type);

    bool SetTrueMode( int &true_mode);
  
    //void WriteOutputFile();
    void WriteTree();
    void CloseOutputFile();
    TFile* GetOutputFile(){return ofile;};
    //void FillBasicData();

    float GetGeneralizedMomentum(){ return  os->amom  ; }
    float GetLeptonZenith()       { return -os->dir[2]; } 

    float GetMGMRE()        { return llMGMRE       ; }
    float GetMGMREnue()     { return llMGMRENUE    ; }
    float GetMGMREnuebar()  { return llMGMRENUEBar ; }
    float GetPi0LL      ()  { return pi0ll         ; }

    float GetProb();
    double Oscillate();
    double Neighbor3D();
    double GetMCweight();
    double GetHondaFluxRatio( int NuType );

    std::string GetMode() { return process_mode ;} 

    int   GetNDecayE()      { return nDecayE; } 

    float   GetLLmpid() { return llmpid   ;  }
    float   GetLLmue () { return llmue    ;  }
    float   GetLLfmom() { return llfmom   ;  } 
    float   GetLLdpos() { return lldpos   ;  }

    void    UseFiTQun( bool x ) { kUseFiTQun = x ; }
    void    UseTauNN( bool x )  { kUseTauNN  = x ; }
    void    DebugMode( bool x )  { kDebugMode  = x ; }
    void    CheckBkg( bool x )  { kCheckBkg  = x ; }
    void    AllHist( bool x )  { kAllHist  = x ; }
    void    MakeNtuple( bool x )  { kMakeNtuple  = x ; }
    void    FermiMotion( bool x ) {kFermiMotion = x; }
    void    OutsideSR( bool x ) {kOutsideSR = x; }
    void    SystNtag( int x ) {kSystNtag = x; }
    void    CorrelatedDecay( int x ) {kCorrelatedDecay = x; }
    void    EnergyScale( int x ) {kEnergyScale = x; }

    void FillNtuple();
    void MakeNtuple();
    
    void BuildEnergyFriend(int type);
    int GetZBin(float zenith, int type);

    bool    SubGeVElikeFQNCRejection();

    //float   GetOscProb();
    ThreeFlavorOscillator *theOscillator;

    fort_func *m_fort;
    BargerPropagator *bNu;

    TTree*   SendTree(std::string treename);
    void   ReadTree(TTree* tree);

    float CalcDrRingVector(int id_ring, int id_vector);
    float CalcDrRingVector2(int id_ring, int id_vector2);
    float CalcDrVector2Vector2(int id1_vector2, int id2_vector2);
    float CalcDrVectorVector(int id1_vector, int id2_vector);
    float CalcOpeningAngle(int id, float mom);//0:electron 1:muon
    float CalcEnergyVector(int id);

    int GetParType(int ring_id);

    int r_max,mu_max,m_max;
    bool pass_cut[11][5];
    int total_with_ntag[3][2][5];
    int lowerSR_without_ntag[3][2][5];
    int higherSR_without_ntag[3][2][5];
    int total;
    float closest_mass_pi0_reco,total_mass,two_elike_mass,total_mom,all_ring_mass,all_ring_mom,all_mulike_mass,all_mulike_mom,total_distance;
    float weight,mc_weight,osc_weight;
    int sample_num,event_num;
    int event_type,nPar,nPar2,nRing,n_elike,n_mulike,interaction_type,nNeutron,true_mode;
    int n_elike_pattern,n_elike_angle,n_mulike_pattern,n_mulike_angle;
    bool is_free_proton;
    int n_free_proton,n_true_neutron,n_tagged_neutron_exp[11][100],n_true_decayE;
    float expected_3ring_events_electron,expected_3ring_events_muon;
    //count p->eee event type
    int n_eee,n_eeeg,n_eeep,n_eeen,n_eeegp,n_eeenp,n_eeegn;
    //count p->epi event type
    int n_noint,n_abs,n_scat,n_charge,n_prod,n_match_e,n_match_1gamma,n_match_2gamma,n_match_0gamma;
    //count p->mupi event type
    int n_match_mu;
    float total_mass_high,total_mass_low,total_mom_high,total_mom_low,pi_mass_low,pi_mass_high;
    
  private:

    float   pi0ll    ;
    float   llmpid   ; 
    float   llmue    ; 
    float   llfmom   ; 
    float   lldpos   ; 

    float llMGMRE       ;
    float llMGMRENUE    ;
    float llMGMRENUEBar ;

    void Initialize();
    void ZeroStructure();
    void Process_pepi();
    void Process_pmupi();
    void Process_peee();
    void Process_pmumumu();
    void Process_pemumu();
    void Process_pmuee();
    void Process_peemu();
    void Process_single();
    void Process_subgev_multiring();
    void Process_subgev_onemulike();
    void Process_subgev_oneelike();
    void Process_cosmic_muon();
    TLorentzVector GetTLorentzVectorVector(int index);
    TLorentzVector GetTLorentzVectorVector2(int index);
    TLorentzVector GetTLorentzVectorRing(int index, int type);
    TVector3 GetTVectorVector(int index);
    TVector3 GetTVectorVector2(int index);
    TVector3 GetTVectorRing(int index);
    TVector3 GetTVectorRingPID(int index, int id);//0 e-like, 1 mu-like

    // event classification functions
    bool SetEventTypeFC(int &);

    DataManager * dm;
    TTree * _ltree;

    // for Multi-GeV Multi-Ring electron likelihood 
    LikelihoodManager * llm;
//  Pi0Likelihood   * llpi0;

    TArrayF * fsiweight;

    // mode related
    std::string process_mode;
    std::string process_input;
    std::string output_ntuple;
    std::string output_hist;
    float live_time_weight;
    unsigned int dataflag;
    unsigned int skgen;


    // related to output
    TFile * ofile;
    TTree * otree;
    TTree * otree2;
    TTree * otree3;
    outputOscStructure * os;
    //values for ntuple
    int o_ipnu;
    int o_itype;
    float o_pnu;
    float o_dir[3];
    float o_amom;
    int o_mode;
    int o_nring;
    int o_nmulike;
    float o_total_mass;
    float o_total_mom;
    float o_dlfct;
    float o_vertex_x;
    float o_vertex_y;
    float o_vertex_z;
    float o_prob_angle[5];
    float o_probms_e[5];
    float o_probms_mu[5];
    float o_prob_pattern[5];
    float o_prmslg_e[5];
    float o_prmslg_mu[5];
    float o_mmom[5];
    float o_dir_x[5];
    float o_dir_y[5];
    float o_dir_z[5];
    float o_ang[5];
    float o_ange[5];
    float o_angm[5];
    float o_mmom_min;
    float o_mmom_mid;
    float o_mmom_max;

  // related to event classification
  int   nDecayE ; 
  float ptot[3];
  unsigned int   nhitac_cut[4];
  unsigned int ehit_cut_1[4];
  unsigned int ehit_cut_2[4];
  float potot_cut[4];
  int graph_point,graph_point_2;
 

  // Private Access to tree elements

  //NTAG
  TypeWrapper<UInt_t> ntag_nn;
  TypeWrapper<Int_t> ntag_mctruth_nn;
  //------------------------------
    
  TypeWrapper<Int_t>  nev;
  TypeWrapper<Int_t>  nsub;
  TypeWrapper<Int_t>  date;
  TypeWrapper<Int_t>  time;
  TypeWrapper<Float_t> pnu;
  TypeWrapper<Float_t> dirnu;
  TypeWrapper<Int_t>   ipnu;
  TypeWrapper<Int_t>   mode;
  TypeWrapper<Float_t> evis;
  TypeWrapper<Float_t> wall;
  TypeWrapper<Float_t> pmomv;
  TypeWrapper<Float_t> pmomv2;
  TypeWrapper<Int_t>   numnu;
  TypeWrapper<Int_t>   nring;
  TypeWrapper<Int_t>   npar;
  TypeWrapper<Int_t>   npar2;
  TypeWrapper<Int_t>   Iflvc;
  TypeWrapper<Int_t>   iprntidx;
  TypeWrapper<Int_t>   Npvc;
  TypeWrapper<Int_t>   Ipvc;
  TypeWrapper<Int_t>   Iorgvc;
  TypeWrapper<Int_t>   Ichvc;
  TypeWrapper<Float_t>   Abspvc;
  TypeWrapper<Float_t>   pscnd;
  TypeWrapper<Int_t>   nscndprt;
  TypeWrapper<Int_t>   iprtscnd;
  TypeWrapper<Float_t>   tscnd;
  TypeWrapper<Int_t>   iprntprt;
  TypeWrapper<Int_t>   lmecscnd;
  TypeWrapper<Int_t>   iprnttrk;
  TypeWrapper<UInt_t>  nhitac;
  TypeWrapper<UInt_t>  nsube;
  TypeWrapper<UInt_t>  ndcy;
  TypeWrapper<UInt_t>  ngate;
  TypeWrapper<UInt_t>  nbye;
  TypeWrapper<Float_t> potot;
  TypeWrapper<UInt_t>  ip; 
  TypeWrapper<UInt_t>  ipv; 
  TypeWrapper<UInt_t>  ipv2; 
  TypeWrapper<UInt_t>  iorg; 
  TypeWrapper<Float_t> pos;
  TypeWrapper<Float_t> posv;
  TypeWrapper<Float_t> ang;
  TypeWrapper<Float_t> ange;
  TypeWrapper<Float_t> angm;
  TypeWrapper<Float_t> amom;
  TypeWrapper<Float_t> amome;
  TypeWrapper<Float_t> amomm;
  TypeWrapper<Float_t> prmslg;
  TypeWrapper<Float_t> probms;
  TypeWrapper<Float_t> Dlfct;
  TypeWrapper<Float_t> dir;
  TypeWrapper<Float_t> dirv;
  TypeWrapper<Float_t> dirv2;
  TypeWrapper<Float_t> dirtotmue;
  TypeWrapper<Float_t> etotmue;
  TypeWrapper<Float_t> msdir;
  TypeWrapper<Int_t> catpc;
  TypeWrapper<Int_t> catpc_qrat_corrected;
  TypeWrapper<Int_t> nmue;
  TypeWrapper<UInt_t> etype;
  TypeWrapper<Float_t> ehit;
  TypeWrapper<Float_t> etime;
  TypeWrapper<Float_t> egood;

  TypeWrapper<UInt_t> hstate;
  TypeWrapper<Float_t> wallv;
  TypeWrapper<Float_t> wallv2;
  TypeWrapper<Float_t>  fit_dir;
  TypeWrapper<Int_t>    fit_pid;
  TypeWrapper<Float_t>  fit_len;
  TypeWrapper<Float_t>  fit_mom;
  TypeWrapper<Int_t>  um_ehit8m;
  TypeWrapper<Int_t>  um_ohit8m;

  TypeWrapper<Int_t>    sh_id  ;
  TypeWrapper<Float_t>  sh_mean;
  TypeWrapper<Float_t>  sh_delta;
  TypeWrapper<Float_t>  sh_meanq;
  TypeWrapper<Float_t>  sh_chi1p;

  TypeWrapper<Float_t> flxb03;
  TypeWrapper<Float_t> flxh06;
  TypeWrapper<Float_t> flxh11;

  TypeWrapper<Float_t> oscwgt;

  // fiTQun related
  TypeWrapper<Float_t> fqpi0nll;
  TypeWrapper<Float_t> fq1rnll;
  TypeWrapper<Float_t> fqpi0mass;
  TypeWrapper<Float_t> fq4rnll;

  bool kUseFiTQun;
  bool kUseTauNN;
  bool kDebugMode;
  bool kCheckBkg;
  bool kAllHist;
  bool kMakeNtuple;
  bool kFermiMotion;
  bool kOutsideSR;
  int kSystNtag;
  int kCorrelatedDecay;
  int kEnergyScale;

  //Tau_NN related
  TypeWrapper<Double_t> NN_output;
  TypeWrapper<Int_t>    NN_selected;

  //friend file related
  TypeWrapper<Float_t> ErmsHax;
  TypeWrapper<Int_t> nEAveHax;

  int Good;
  int Bad;


};


#endif


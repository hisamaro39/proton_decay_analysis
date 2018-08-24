#include "TH2F.h"
#include "TTree.h"
#include "OscNtupleManager.h"
#include "outputOscStructure.h"
#include "osc_types.h"

const Int_t nbins_rec_type=4;
const Int_t nbins_nn=30;
const Int_t nbins_muedk=10;
const Int_t nbins_type=EndOfTypes-1;
const Float_t Elow_resonance=2.5;
const Float_t Ehigh_resonance=11.22;

class Energy_recon{
 private:
  std::string ref_file;
  TFile* ofile;
  TTree* osc_tuple_origin;
  TTree* osc_tuple_rec;
  TTree* otree_erecon;
  TH1F* amom_offset_hist[nbins_rec_type][nbins_nn][nbins_muedk];
  Float_t fit_params[nbins_rec_type][nbins_nn][nbins_muedk][5];
  Float_t amom_offset[nbins_rec_type][nbins_nn][nbins_muedk];
  outputOscStructure * os2;
  void SetOutput_ERecon();
  void SetBranch_ERecon();
  
 public:
  Energy_recon();
  ~Energy_recon();
  void ReadTree(OscNtupleManager* om);
  void SendTree(OscNtupleManager* om);
  void DeleteTree();
  void make_table();
  void remake_tree();
  void offset_fix(Bool_t read_io);
  void setFileName(std::string file){ref_file=file;};
  std::string getFileName(){return ref_file;};
  Double_t MAD(TH1F* h1, Double_t median_value);
  Int_t rec_type[nbins_type];
  void GetSeparation(std::string cardname);
  void SaveFunction(std::string func_file);
  void GetFunction(std::string func_file);
  void SetOutputFile(TFile* file);
  
};

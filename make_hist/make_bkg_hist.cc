#include <vector>
void make_bkg_hist(){
  string type[] = {"p_epi","p_mupi","p_eee","p_mumumu","p_emumu","p_muee","fcmc","fcdt"};
  int sample_mode = 2;
  int neut_mode = 0;//0:nue 1:numu
  int int_mode = 1;//interaction mode. Negative means anti neutrino.

  int nring=-1,mulike=-1,michel=-1;
  if(sample_mode==2){
    nring=1;
    mulike=0;
    michel=0;
  }
  if(sample_mode==5){
    nring=1;
    mulike=0;
    michel=1;
  }

  gROOT->SetStyle("Plain");
  gStyle->SetOptStat(0);
  gStyle->SetPalette(1);
  vector<string> hist;
  vector<int> rebin,logy;

  //hist list
  /*hist.push_back("nRing_cut1");  
  rebin.push_back(1); logy.push_back(0);
  hist.push_back("nElikeRing_angle_nring3_cut1");  
  rebin.push_back(1);logy.push_back(0);
  hist.push_back("nElikeRing_angle_nring2_cut1");  
  rebin.push_back(1);logy.push_back(0);
  hist.push_back("nMulikeRing_angle_nring3_cut1");  
  rebin.push_back(1);logy.push_back(0);
  hist.push_back("nMulikeRing_angle_nring2_cut1");  
  rebin.push_back(1);logy.push_back(0);
  hist.push_back("n_michel_electron_cut3");  
  rebin.push_back(1);logy.push_back(0);
  hist.push_back("ntag_multiplicity_cut4");  
  rebin.push_back(1);logy.push_back(0);
  hist.push_back("mass_proton_reco_cut4");  
  rebin.push_back(5);logy.push_back(0);
  hist.push_back("mom_proton_reco_cut4");  
  rebin.push_back(5);logy.push_back(0);*/
  hist.push_back("mass_two_elike_reco_cut4");  
  rebin.push_back(5);logy.push_back(0);

  TFile *input;
  TH1 *this_hist;
  for(int h=0;h<hist.size();h++){
    input = TFile::Open(Form("../output/fcmc.sk4.mode_%s_pi0cut_check_bkg.root",type[sample_mode].c_str()));
    TCanvas *c = new TCanvas(Form("canvas%d",h),"",800,600);
    if(logy[h]) c->SetLogy();
    string save_name = "";
    if(int_mode>0) this_hist = (TH1*) input->Get(Form("%s_nring%d_mulike%d_michel%d_type%d_mode_pos%d",hist[h].c_str(),nring,mulike,michel,neut_mode,int_mode));
    else this_hist = (TH1*) input->Get(Form("%s_nring%d_mulike%d_michel%d_type%d_mode_neg%d",hist[h].c_str(),nring,mulike,michel,neut_mode,abs(int_mode)));
    this_hist->Rebin(rebin[h]);
    this_hist->SetLineWidth(2);
    this_hist->Draw("hist E0");
    if(int_mode>0) c->SaveAs(Form("hist_bkg/single_%s_nring%d_mulike%d_michel%d_input_fcmc_mode_%s_type%d_mode_pos%d.pdf",hist[h].c_str(),nring,mulike,michel,type[sample_mode[h]].c_str(),neut_mode,int_mode));
    else c->SaveAs(Form("hist_bkg/single_%s_nring%d_mulike%d_michel%d_input_fcmc_mode_%s_type%d_mode_neg%d.pdf",hist[h].c_str(),nring,mulike,michel,type[sample_mode[h]].c_str(),neut_mode,abs(int_mode)));
  }

}

#include <vector>
void make_single_plot_th1(){
  string type[] = {"p_epi","p_mupi","p_eee","p_mumumu","p_emumu","p_muee","fcmc","fcdt"};

  gROOT->SetStyle("Plain");
  gStyle->SetOptStat(0);
  gStyle->SetPalette(1);
  vector<string> hist;
  vector<int> input_type,mode_type;
  hist.clear();input_type.clear();mode_type.clear();
  //hist list
  hist.push_back("nRing_cut1_nring0_mulike0_michel0");  
  input_type.push_back(5);mode_type.push_back(5);
  hist.push_back("nMulikeRing_nring3_cut1_nring0_mulike0_michel0");  
  input_type.push_back(5);mode_type.push_back(5);
  hist.push_back("nMulikeRing_nring2_cut1_nring0_mulike0_michel0");  
  input_type.push_back(5);mode_type.push_back(5);

  TFile *input;
  for(int h=0;h<hist.size();h++){
    input = TFile::Open(Form("../output/%s.sk4.mode_%s.root",type[input_type[h]].c_str(),type[mode_type[h]].c_str()));
    TCanvas *c = new TCanvas(Form("canvas%d",h),"",800,600);
    string save_name = "";
    TH2* this_hist = (TH2*) input->Get(hist[h].c_str());
    this_hist->SetLineWidth(2);
    this_hist->Draw("hist");
    c->SaveAs(Form("hist/single_%s_input_%s_mode_%s.pdf",hist[h].c_str(),type[input_type[h]].c_str(),type[mode_type[h]].c_str()));
  }

}

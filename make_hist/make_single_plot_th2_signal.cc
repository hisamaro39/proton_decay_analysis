#include <vector>
void make_single_plot_th2_signal(){
  string type[] = {"p_epi","p_mupi","p_eee","p_mumumu","p_emumu","p_muee","fcmc","fcdt"};
  int mode_id = 1;

  gROOT->SetStyle("Plain");
  gStyle->SetOptStat(0);
  gStyle->SetPalette(1);
  vector<string> hist;
  vector<int> input_type,mode_type;
  hist.clear();input_type.clear();mode_type.clear();
  //hist list
  hist.push_back("mass_mom_proton_reco_cut6_nring0_mulike0_michel0");  
  input_type.push_back(mode_id);mode_type.push_back(mode_id);

  TBox *box = new TBox(800,0,1050,100);
  box->SetFillStyle(0);
  box->SetLineWidth(2);
  box->SetLineColor(1);
  TBox *box2 = new TBox(800,100,1050,250);
  box2->SetFillStyle(0);
  box2->SetLineWidth(2);
  box2->SetLineColor(1);
  TFile *input;
  for(int h=0;h<hist.size();h++){
    input = TFile::Open(Form("../output/%s.sk4.mode_%s.root",type[input_type[h]].c_str(),type[mode_type[h]].c_str()));
    TCanvas *c = new TCanvas(Form("canvas%d",h),"",800,600);
    string save_name = "";
    TH2* this_hist_fp0 = (TH2*) input->Get(Form("%s_fp0",hist[h].c_str()));
    TH2* this_hist_fp1 = (TH2*) input->Get(Form("%s_fp1",hist[h].c_str()));
    this_hist_fp0->SetMinimum(0);
    this_hist_fp1->SetMinimum(0);
    //this_hist_fp0->SetMarkerStyle(8);
    //this_hist_fp1->SetMarkerStyle(8);
    this_hist_fp1->SetMarkerColor(2);
    this_hist_fp0->Draw();
    this_hist_fp1->Draw("same");
    box->Draw();
    box2->Draw();
    c->SaveAs(Form("hist/single_%s_input_%s_mode_%s.pdf",hist[h].c_str(),type[mode_id].c_str(),type[mode_id].c_str()));
  }

}

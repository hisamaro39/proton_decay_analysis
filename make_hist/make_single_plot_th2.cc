#include <vector>
void make_single_plot_th2(){
  string type[] = {"p_epi","p_mupi","p_eee","p_mumumu","p_emumu","p_muee","fcmc","fcdt"};

  gROOT->SetStyle("Plain");
  gStyle->SetOptStat(0);
  gStyle->SetPalette(1);
  vector<string> hist;
  vector<int> input_type,mode_type;
  hist.clear();input_type.clear();mode_type.clear();
  //hist list
  hist.push_back("mass_mom_proton_reco_weight_osc_cut4_nring0_mulike1_michel1");  
  input_type.push_back(6);mode_type.push_back(5);
  hist.push_back("mass_mom_proton_reco_weight_osc_cut4_nring1_mulike1_michel1");  
  input_type.push_back(6);mode_type.push_back(5);

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
    TH2* this_hist = (TH2*) input->Get(hist[h].c_str());
    this_hist->SetMinimum(0);
    this_hist->Draw("colz");
    box->Draw();
    box2->Draw();
    c->SaveAs(Form("hist/single_%s_input_%s_mode_%s.pdf",hist[h].c_str(),type[input_type[h]].c_str(),type[mode_type[h]].c_str()));
  }

}

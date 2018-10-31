#include <vector>
#include "TGraphAsymmErrors.h"
void compare_library(){
  string type[] = {"p_epi","p_mupi","p_eee","p_mumumu","p_emumu","p_muee","p_eemu","fcmc","fcdt"};

  gROOT->SetStyle("Plain");
  gStyle->SetOptStat(0);
  vector<string> hist;
  vector<vector<string>> hist_set;
  vector<vector<int>> input_type_set,mode_type_set;
  vector<int> scale,dology,input_type,add_ratio,mode_type,use_validation;
  hist.clear();hist_set.clear();scale.clear();dology.clear();input_type.clear();
  input_type_set.clear();add_ratio.clear();mode_type_set.clear();mode_type.clear();use_validation.clear();

  hist.push_back("nRing_cut1_nring0_mulike0_michel0");
  input_type.push_back(7);
  mode_type.push_back(0);
  dology.push_back(0);

  hist.push_back("cut_flow_nring2_mulike0_michel0");
  input_type.push_back(7);
  mode_type.push_back(0);
  dology.push_back(0);

  TFile *input_old,*input_new;
  TH1 *first_hist;
  for(int h=0;h<hist.size();h++){
    TCanvas *c = new TCanvas(Form("canvas%d",h),"",800,600);
    if(dology[h]) c->SetLogy();
    string save_name = "";
    save_name = "hist/compare_library_";
    input_old = TFile::Open(Form("../output/%s.sk4.mode_%s.root",type[input_type[h]].c_str(),type[mode_type[h]].c_str()));
    input_new = TFile::Open(Form("/home/matanaka/disk02_usr6/proton_decay/proton_decay_analysis/output/%s.sk4.mode_%s.root",type[input_type[h]].c_str(),type[mode_type[h]].c_str()));
    TH1* hist_old = (TH1*) input_old->Get(hist[h].c_str());
    TH1* hist_new = (TH1*) input_new->Get(hist[h].c_str());
    hist_old->SetLineWidth(2);
    hist_old->SetLineColor(1);
    hist_new->SetLineWidth(2);
    hist_new->SetLineColor(2);
    hist_old->Draw("hist E0");
    hist_new->Draw("hist E0 same");
    save_name += hist[h];
    save_name += "_" + type[input_type[h]] + "_" + type[mode_type[h]];
    save_name += ".pdf";
    cout << save_name << endl;
    c->SaveAs(save_name.c_str());
    //c->SaveAs("hist/temp.pdf");

  }

}

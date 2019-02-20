#include <vector>
#include "TGraphAsymmErrors.h"
void make_compare_plot_th2(){
  string type[] = {"p_epi","p_mupi","p_eee","p_mumumu","p_emumu","p_muee","p_eemu","fcmc","fcdt"};

  gROOT->SetStyle("Plain");
  gStyle->SetOptStat(0);
  vector<string> hist;
  vector<vector<string>> hist_set;
  vector<vector<int>> input_type_set,mode_type_set;
  vector<int> scale,dology,input_type,add_ratio,mode_type,use_validation;
  hist.clear();hist_set.clear();scale.clear();dology.clear();input_type.clear();
  input_type_set.clear();add_ratio.clear();mode_type_set.clear();mode_type.clear();use_validation.clear();

  hist.push_back("residual_total_mass_vs_true_min_mom_muon_nring3_mulike1_fp1");
  hist.push_back("residual_total_mass_vs_true_min_mom_muon_nring3_mulike2_fp1");
  hist.push_back("residual_total_mass_vs_true_min_mom_muon_nring3_mulike3_fp1");
  input_type.push_back(5);
  input_type.push_back(4);
  input_type.push_back(3);
  mode_type.push_back(5);
  mode_type.push_back(4);
  mode_type.push_back(3);
  add_ratio.push_back(0);scale.push_back(1);use_validation.push_back(1);
  hist_set.push_back(hist);dology.push_back(0);hist.clear();input_type_set.push_back(input_type);input_type.clear();
  mode_type_set.push_back(mode_type);mode_type.clear();

  TFile *input;
  for(int s=0;s<hist_set.size();s++){
    string save_name = "";
    save_name = "hist/compare_";
    for(int h=0;h<hist_set[s].size();h++){
      if(use_validation[s]) input = TFile::Open(Form("../output/%s.sk4.mode_%s_validation.root",type[input_type_set[s].at(h)].c_str(),type[mode_type_set[s].at(h)].c_str()));
      else input = TFile::Open(Form("../output/%s.sk4.mode_%s.root",type[input_type_set[s].at(h)].c_str(),type[mode_type_set[s].at(h)].c_str()));
      TH2* this_hist = (TH2*) input->Get(hist_set[s].at(h).c_str());
      if(h==0) this_hist->Draw("colz");
      else this_hist->Draw("same colz");
      save_name += hist_set[s].at(h);
      save_name += "_" + type[input_type_set[s].at(h)] + "_" + type[input_type_set[s].at(h)];
    }
    save_name += ".pdf";
    //c->SaveAs(save_name.c_str());
    //c->SaveAs("hist/temp.pdf");

  }

}

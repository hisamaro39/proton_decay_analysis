#include <vector>
void make_single_plot_th1_apfit(){

  gROOT->SetStyle("Plain");
  gStyle->SetOptStat(0);
  gStyle->SetPalette(1);
  vector<string> hist;
  vector<int> input_type,mode_type,rebin;
  hist.clear();input_type.clear();mode_type.clear();rebin.clear();
  //hist list
  /*hist.push_back("h_diff_dlfct");  
  hist.push_back("h_diff_total_mass");  
  hist.push_back("h_diff_total_mom");  
  hist.push_back("h_diff_mmom");  
  hist.push_back("h_diff_prob_angle");  
  hist.push_back("h_diff_mmom_min");  
  hist.push_back("h_diff_mmom_mid");  
  hist.push_back("h_diff_mmom_max");*/  
  //hist.push_back("h_diff_vertex_x");  
  //hist.push_back("h_diff_vertex_y");  
  //hist.push_back("h_diff_vertex_z");  
  //hist.push_back("h_diff_dlfct_same_vertex");  
  //hist.push_back("h_diff_dlfct_different_vertex");  
  //hist.push_back("h_diff_prob_angle_same_vertex");  
  //hist.push_back("h_diff_prob_angle_different_vertex");  
  //hist.push_back("h_diff_probms_e");  
  //hist.push_back("h_diff_probms_mu");  
  hist.push_back("h_diff_mmom_same_vertex");  
  hist.push_back("h_diff_mmom_different_vertex");  

  TFile *input = TFile::Open("output/output_compare_apfit_p_mumumu.root");
  for(int h=0;h<hist.size();h++){
    TCanvas *c = new TCanvas(Form("canvas%d",h),"",800,600);
    c->SetLogy();
    string save_name = "";
    TH1* this_hist = (TH1*) input->Get(hist[h].c_str());
    this_hist->SetMinimum(1);
    this_hist->SetLineWidth(2);
    this_hist->Draw("hist E0");
    c->SaveAs(Form("hist/compare_apfit_p_mumumu_%s.pdf",hist[h].c_str()));
  }

}

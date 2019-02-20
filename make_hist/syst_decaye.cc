#include <vector>
void syst_decaye(){
  string mode = "onemulike";

  gROOT->SetStyle("Plain");
  gStyle->SetOptStat(0);

  cout << "mode=" << mode << endl;
  TH1 *first_hist;
  TFile *input_mc = TFile::Open(Form("../output/fcmc.sk4.mode_subgev_%s.root",mode.c_str()));//mc
  TFile *input_data = TFile::Open(Form("../output/fcdt.sk4.mode_subgev_%s.root",mode.c_str()));//mc
  TH1* hist_mc = (TH1*) input_mc->Get("n_michel_electron");
  TH1* hist_data = (TH1*) input_data->Get("n_michel_electron");
  hist_data->SetLineWidth(2);
  hist_mc->SetLineWidth(2);
  cout << "entries data/mc=" << hist_data->GetEntries() << "/" << hist_mc->GetEntries() << endl;
  TCanvas *c1 = new TCanvas("c1","",800,600);
  hist_data->Draw();
  hist_mc->SetLineColor(2);
  hist_mc->Draw("same hist");
  c1->SaveAs(Form("hist/compare_data_mc_n_michel_electron_%s.pdf",mode.c_str()));
  float evt_data = hist_data->Integral();
  float evt_mc = hist_mc->Integral();
  float norm_factor = evt_data/evt_mc;
  cout << "total events data/mc=" << evt_data << "/" << evt_mc << endl;
  cout << "norm_factor=" << norm_factor << endl;
  float zero_decaye_data = hist_data->GetBinContent(1);
  float zero_decaye_mc = hist_mc->GetBinContent(1);
  float one_decaye_data = hist_data->GetBinContent(2);
  float one_decaye_mc = hist_mc->GetBinContent(2);
  if(mode=="onemulike") cout << "one decayE ecvents data/mc=" << one_decaye_data << "/" << one_decaye_mc << endl;
  if(mode=="oneelike") cout << "zero decayE ecvents data/mc=" << zero_decaye_data << "/" << zero_decaye_mc << endl;
  TH1* hist_mc_norm = (TH1*) hist_mc->Clone("hist_mc_norm");
  hist_mc_norm->Scale(norm_factor);
  hist_mc_norm->SetLineColor(4);
  hist_mc_norm->SetLineWidth(2);
  TCanvas *c2 = new TCanvas("c2","",800,600);
  hist_data->Draw();
  hist_mc_norm->Draw("same hist");
  c2->SaveAs(Form("hist/compare_data_mc_n_michel_electron_scaled_%s.pdf",mode.c_str()));
  float evt_mc_norm = hist_mc_norm->Integral();
  float zero_decaye_mc_norm = hist_mc_norm->GetBinContent(1);
  float one_decaye_mc_norm = hist_mc_norm->GetBinContent(2);
  float err_one_decaye_mc_norm = hist_mc_norm->GetBinError(2);
  float err_zero_decaye_mc_norm = hist_mc_norm->GetBinError(1);
  cout << "total events data/scale_mc=" << evt_data << "/" << evt_mc_norm << endl;
  if(mode=="onemulike"){
    cout << "one decayE ecvents data=" << one_decaye_data << endl;
    cout << "one decayE ecvents norm_mc=" << one_decaye_mc_norm << " +- " << err_one_decaye_mc_norm << endl;
    float diff = one_decaye_data - one_decaye_mc_norm;
    float diff_err = err_one_decaye_mc_norm;
    float diff_ratio = diff / one_decaye_data;
    float diff_ratio_err = err_one_decaye_mc_norm / one_decaye_data; 
  }
  if(mode=="oneelike"){
    cout << "zero decayE ecvents data=" << zero_decaye_data << endl;
    cout << "zero decayE ecvents norm_mc=" << zero_decaye_mc_norm << " +- " << err_zero_decaye_mc_norm << endl;
    float diff = zero_decaye_data - zero_decaye_mc_norm;
    float diff_err = err_zero_decaye_mc_norm;
    float diff_ratio = diff / zero_decaye_data;
    float diff_ratio_err = err_zero_decaye_mc_norm / zero_decaye_data; 
  }
  cout << "data - mc = " << diff << " +- " << diff_err << endl;
  cout << "diff / data = " << diff_ratio << " +- " << diff_ratio_err << endl; 

}

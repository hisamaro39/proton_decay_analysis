#include <vector>
void syst_pid3(){

  gROOT->SetStyle("Plain");
  gStyle->SetOptStat(0);

  TH1 *first_hist;
  TFile *input_mc = TFile::Open("../output/fcmc.sk4.mode_subgev_multiring.root");//mc
  TFile *input_data = TFile::Open("../output/fcdt.sk4.mode_subgev_multiring.root");//data
  TH1* hist_mc = (TH1*) input_mc->Get("nElikeRing_angle_nring3");
  TH1* hist_data = (TH1*) input_data->Get("nElikeRing_angle_nring3");
  hist_data->SetLineWidth(2);
  hist_mc->SetLineWidth(2);
  cout << "entries data/mc=" << hist_data->GetEntries() << "/" << hist_mc->GetEntries() << endl;
  TCanvas *c1 = new TCanvas("c1","",800,600);
  hist_data->Draw();
  hist_mc->SetLineColor(2);
  hist_mc->Draw("same hist");
  c1->SaveAs("hist/compare_data_mc_nElikeRing_angle_nring3_subgev_multiring.pdf");
  float evt_data = hist_data->Integral();
  float evt_mc = hist_mc->Integral();
  float norm_factor = evt_data/evt_mc;
  cout << "total events data/mc=" << evt_data << "/" << evt_mc << endl;
  cout << "norm_factor=" << norm_factor << endl;
  float zero_elike_data = hist_data->GetBinContent(1);
  float one_elike_data = hist_data->GetBinContent(2);
  float two_elike_data = hist_data->GetBinContent(3);
  float three_elike_data = hist_data->GetBinContent(4);
  TH1* hist_mc_norm = (TH1*) hist_mc->Clone("hist_mc_norm");
  hist_mc_norm->Scale(norm_factor);
  hist_mc_norm->SetLineColor(4);
  hist_mc_norm->SetLineWidth(2);
  TCanvas *c2 = new TCanvas("c2","",800,600);
  hist_data->Draw();
  hist_mc_norm->Draw("same hist");
  c2->SaveAs("hist/compare_data_mc_nElikeRing_angle_nring3_subgev_multiring_scaled.pdf");
  float zero_elike_mc = hist_mc_norm->GetBinContent(1);
  float one_elike_mc = hist_mc_norm->GetBinContent(2);
  float two_elike_mc = hist_mc_norm->GetBinContent(3);
  float three_elike_mc = hist_mc_norm->GetBinContent(4);
  cout << "zero e-like events data/mc =" << zero_elike_data << "/" << zero_elike_mc << endl;
  cout << "one e-like events data/mc =" << one_elike_data << "/" << one_elike_mc << endl;
  cout << "two e-like events data/mc =" << two_elike_data << "/" << two_elike_mc << endl;
  cout << "three e-like events data/mc =" << three_elike_data << "/" << three_elike_mc << endl;
  float err_zero_elike = sqrt(zero_elike_data)/zero_elike_data; 
  float err_one_elike = sqrt(one_elike_data)/one_elike_data; 
  float err_two_elike = sqrt(two_elike_data)/two_elike_data; 
  float err_three_elike = sqrt(three_elike_data)/three_elike_data; 
  float diff_zero_elike = (zero_elike_data - zero_elike_mc) / zero_elike_data;
  float diff_one_elike = (one_elike_data - one_elike_mc) / one_elike_data;
  float diff_two_elike = (two_elike_data - two_elike_mc) / two_elike_data;
  float diff_three_elike = (three_elike_data - three_elike_mc) / three_elike_data;

  cout << "Diff zero elike  = " << diff_zero_elike << " +- " << diff_zero_elike * err_zero_elike << endl; 
  cout << "Diff one elike  = " << diff_one_elike << " +- " << diff_one_elike * err_one_elike << endl; 
  cout << "Diff two elike  = " << diff_two_elike << " +- " << diff_two_elike * err_two_elike << endl; 
  cout << "Diff three elike  = " << diff_three_elike << " +- " << diff_three_elike * err_three_elike << endl; 

}

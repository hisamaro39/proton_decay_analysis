#include <vector>
void syst_ring2(){

  gROOT->SetStyle("Plain");
  gStyle->SetOptStat(0);

  TH1 *first_hist;
  TFile *input_mc = TFile::Open("../output/fcmc.sk4.mode_subgev_multiring.root");//mc
  TFile *input_data = TFile::Open("../output/fcdt.sk4.mode_subgev_multiring.root");//data
  TH1* hist_mc = (TH1*) input_mc->Get("nRing");
  TH1* hist_data = (TH1*) input_data->Get("nRing");
  hist_data->SetLineWidth(2);
  hist_mc->SetLineWidth(2);
  cout << "entries data/mc=" << hist_data->GetEntries() << "/" << hist_mc->GetEntries() << endl;
  TCanvas *c1 = new TCanvas("c1","",800,600);
  hist_data->Draw();
  hist_mc->SetLineColor(2);
  hist_mc->Draw("same hist");
  c1->SaveAs("hist/compare_data_mc_nRing_subgev.pdf");
  float evt_data = hist_data->Integral();
  float evt_mc = hist_mc->Integral();
  float norm_factor = evt_data/evt_mc;
  cout << "total events data/mc=" << evt_data << "/" << evt_mc << endl;
  cout << "norm_factor=" << norm_factor << endl;
  float three_ring_data = hist_data->GetBinContent(4);
  TH1* hist_mc_norm = (TH1*) hist_mc->Clone("hist_mc_norm");
  hist_mc_norm->Scale(norm_factor);
  hist_mc_norm->SetLineColor(4);
  hist_mc_norm->SetLineWidth(2);
  TCanvas *c2 = new TCanvas("c2","",800,600);
  hist_data->Draw();
  hist_mc_norm->Draw("same hist");
  c2->SaveAs("hist/compare_data_mc_nRing_subgev_scaled.pdf");
  float three_ring_mc = hist_mc_norm->GetBinContent(4);
  cout << "three ring events data/mc =" << three_ring_data << "/" << three_ring_mc << endl;
  float err_three_ring = sqrt(three_ring_data)/three_ring_data; 
  float diff_three_ring = (three_ring_data - three_ring_mc) / three_ring_data;

  cout << "Diff three ring  = " << diff_three_ring << " +- " << diff_three_ring * err_three_ring << endl; 

}

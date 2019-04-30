#include <vector>
void syst_decaye_2(){
  int outsideSR = 1;
  //string mode = "subgev_multiring";
  string mode = "p_mumumu";
  int period = 4;//5:sk1-3
  //string hist_name = "n_michel_electron_nring3_mulike3";
  string hist_name = "n_michel_electron_cut3_nring0_mulike0_michel3";
  int michel = 3;//3:2or3

  gROOT->SetStyle("Plain");
  gStyle->SetOptStat(0);

  cout << "mode=" << mode << endl;
  cout << "sk" << period << endl;
  TH1 *first_hist;
  TFile *input_mc,*input_data;
  if(outsideSR){
    if(period==5){
      input_mc = TFile::Open(Form("../output/fcmc_final.sk1_3.mode_%s_outsideSR_wo_livetime.root",mode.c_str()));//mc
      input_data = TFile::Open(Form("../output/fcdt_final.sk1_3.mode_%s_outsideSR_wo_livetime.root",mode.c_str()));//mc
    }else{
      input_mc = TFile::Open(Form("../output/fcmc_final.sk%d.mode_%s_outsideSR.root",period,mode.c_str()));//mc
      input_data = TFile::Open(Form("../output/fcdt_final.sk%d.mode_%s_outsideSR.root",period,mode.c_str()));//mc
    }
  }
  else{
    input_mc = TFile::Open(Form("../output/fcmc_final.sk%d.mode_%s.root",period,mode.c_str()));//mc
    input_data = TFile::Open(Form("../output/fcdt_final.sk%d.mode_%s.root",period,mode.c_str()));//mc
  }
  TH1* hist_mc = (TH1*) input_mc->Get(hist_name.c_str());
  TH1* hist_data = (TH1*) input_data->Get(hist_name.c_str());
  hist_data->SetLineWidth(2);
  hist_mc->SetLineWidth(2);
  cout << "entries data/mc=" << hist_data->GetEntries() << "/" << hist_mc->GetEntries() << endl;
  TCanvas *c1 = new TCanvas("c1","",800,600);
  hist_data->Draw();
  hist_mc->SetLineColor(2);
  hist_mc->Draw("same hist");
  //c1->SaveAs(Form("hist/compare_data_mc_%s_%s_sk%d.pdf",hist_name.c_str(),mode.c_str(),period));
  float evt_data = hist_data->Integral();
  float evt_mc = hist_mc->Integral();
  float norm_factor = evt_data/evt_mc;
  cout << "total events data/mc=" << evt_data << "/" << evt_mc << endl;
  cout << "norm_factor=" << norm_factor << endl;
  float decaye_data = hist_data->GetBinContent(michel+1);
  if(michel==3) decaye_data += hist_data->GetBinContent(michel); 
  float decaye_mc = hist_mc->GetBinContent(michel+1);
  if(michel==3) decaye_mc += hist_mc->GetBinContent(michel);
  cout << michel << " decayE ecvents data/mc=" << decaye_data << "/" << decaye_mc << endl;
  TH1* hist_mc_norm = (TH1*) hist_mc->Clone("hist_mc_norm");
  hist_mc_norm->Scale(norm_factor);
  hist_mc_norm->SetLineColor(4);
  hist_mc_norm->SetLineWidth(2);
  //TCanvas *c2 = new TCanvas("c2","",800,600);
  //hist_data->Draw();
  hist_mc_norm->Draw("same hist");
  c1->SaveAs(Form("hist/compare_data_mc_%s_%s_sk%d.pdf",hist_name.c_str(),mode.c_str(),period));
  float evt_mc_norm = hist_mc_norm->Integral();
  float decaye_mc_norm = hist_mc_norm->GetBinContent(michel+1);
  if(michel==3) decaye_mc_norm += hist_mc_norm->GetBinContent(michel);
  float err_decaye_mc_norm = hist_mc_norm->GetBinError(michel+1);
  if(michel==3) sqrt(pow(err_decaye_mc_norm,2) + pow(hist_mc_norm->GetBinContent(michel),2));
  float stat_err = sqrt(decaye_data)/decaye_data;
  cout << "total events data/scale_mc=" << evt_data << "/" << evt_mc_norm << endl;
  cout << michel << " decayE ecvents data=" << decaye_data << " +- " << sqrt(decaye_data) << endl;
  cout << michel << " decayE ecvents norm_mc=" << decaye_mc_norm << " +- " << err_decaye_mc_norm << endl;
  float diff = decaye_data - decaye_mc_norm;
  float diff_err = diff*stat_err;
  float diff_ratio = diff / decaye_data;
  float diff_ratio_err = diff_ratio*stat_err; 
  cout << "data - mc = " << diff << " +- " << diff_err << endl;
  cout << "diff / data = " << diff_ratio << " +- " << diff_ratio_err << endl; 

}

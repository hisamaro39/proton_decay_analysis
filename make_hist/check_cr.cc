#include <vector>
void check_cr(){
  string mode = "p_eee";
  int nring=1;
  int mulike=0;
  int michel=0;
  int period = 5;//5:sk1-4
  float mc_error = 0.3;//systematic error for background

  gROOT->SetStyle("Plain");
  gStyle->SetOptStat(0);
  vector<string> hist;
  vector<vector<string>> hist_set;
  vector<int> dology,dorebin;
  hist.clear();hist_set.clear();dology.clear();dorebin.clear();
  vector<string> hist_name;
  hist_name.clear();

  hist_name.push_back(Form("nRing_cut1_nring%d_mulike%d_michel%d",nring,mulike,michel));
  dology.push_back(0);dorebin.push_back(1);

  hist_name.push_back(Form("nElikeRing_angle_nring3_cut2_nring%d_mulike%d_michel%d",nring,mulike,michel));
  dology.push_back(0);dorebin.push_back(1);

  hist_name.push_back(Form("nMulikeRing_angle_nring3_cut2_nring%d_mulike%d_michel%d",nring,mulike,michel));
  dology.push_back(0);dorebin.push_back(1);

  hist_name.push_back(Form("n_michel_electron_cut3_nring%d_mulike%d_michel%d",nring,mulike,michel));
  dology.push_back(0);dorebin.push_back(1);

  hist_name.push_back(Form("mass_pi0_reco_cut4_nring%d_mulike%d_michel%d",nring,mulike,michel));
  dology.push_back(0);dorebin.push_back(5);

  hist_name.push_back(Form("ntag_multiplicity_cut5_nring%d_mulike%d_michel%d",nring,mulike,michel));
  dology.push_back(0);dorebin.push_back(1);

  hist_name.push_back(Form("mass_proton_reco_cut5_nring%d_mulike%d_michel%d",nring,mulike,michel));
  dology.push_back(0);dorebin.push_back(5);

  hist_name.push_back(Form("mom_proton_reco_cut5_nring%d_mulike%d_michel%d",nring,mulike,michel));
  dology.push_back(0);dorebin.push_back(5);

  TH1 *first_hist;
  TFile *input_data,*input_mc;
  if(period==5){
    input_data = TFile::Open(Form("../output/fcdt_final.sk1_4.mode_%s_CR.root",mode.c_str()));//data
    input_mc = TFile::Open(Form("../output/fcmc_real.sk1_4.mode_%s_CR.root",mode.c_str()));//mc
  }
  else{
    input_data = TFile::Open(Form("../output/fcdt_final.sk%d.mode_%s_CR.root",period,mode.c_str()));//data
    input_mc = TFile::Open(Form("../output/fcmc_real.sk%d.mode_%s_CR.root",period,mode.c_str()));//mc
  }
  for(int s=0;s<hist_name.size();s++){
    TCanvas *c = new TCanvas(Form("canvas%d",s),"",800,600);
    TH1* hist_data = (TH1*) input_data->Get(hist_name[s].c_str());
    TH1* hist_mc = (TH1*) input_mc->Get(hist_name[s].c_str());
    hist_mc->Rebin(dorebin[s]);
    TH1 *hist_mc_error = (TH1*) hist_mc->Clone("hist_mc_error");
    for(int b=0;b<hist_mc_error->GetNbinsX();b++){
      float error = hist_mc_error->GetBinContent(b+1)*mc_error;
      hist_mc_error->SetBinError(b+1,error);
    }
    hist_data->Rebin(dorebin[s]);
    int max = hist_data->GetMaximum()*1.5;
    hist_mc->SetMaximum(max);
    hist_mc->SetLineWidth(2);
    hist_data->SetLineWidth(2);
    hist_mc_error->SetFillColor(1);
    hist_mc_error->SetFillStyle(3003);
    hist_mc->Draw("hist");
    //hist_mc_error->Draw("same E2");
    hist_data->Draw("same E0");
    c->SaveAs(Form("hist/check_cr_%s_%s_sk%d.pdf",hist_name[s].c_str(),mode.c_str(),period));
  }

}

#include "TGraphAsymmErrors.h"
void make_compare_plot_th1_syst_single(){
  string period = "sk4";
  string type[] = {"p_epi","p_mupi","p_eee_final","p_eee_miura","p_mumumu_def","p_mumumu_miura","p_emumu_def","p_emumu_miura","p_muee_def","p_muee_miura","p_eemu","fcmc","fcdt","single_mu_def","single_mu_miura"};
  string syst = "fermigas";

  gROOT->SetStyle("Plain");
  gStyle->SetOptStat(0);
  vector<string> hist;
  vector<int> dology,input_type,add_ratio,mode_type,use_validation,scale,rebin;
  hist.clear();dology.clear();input_type.clear();scale.clear();rebin.clear();
  add_ratio.clear();mode_type.clear();use_validation.clear();

  /*hist.push_back("total_true_mass_fp0");
  input_type.push_back(3);mode_type.push_back(3);scale.push_back(1);
  use_validation.push_back(1);dology.push_back(1);add_ratio.push_back(1);rebin.push_back(1);

  hist.push_back("total_true_mom_fp0");
  input_type.push_back(3);mode_type.push_back(3);scale.push_back(1);
  use_validation.push_back(1);dology.push_back(1);add_ratio.push_back(1);rebin.push_back(1);

  hist.push_back("true_mom_lepton_fp0");
  input_type.push_back(3);mode_type.push_back(3);scale.push_back(1);
  use_validation.push_back(1);dology.push_back(1);add_ratio.push_back(1);rebin.push_back(1);

  hist.push_back("nRing_cut1_nring1_mulike0_michel0");
  input_type.push_back(3);mode_type.push_back(3);scale.push_back(1);
  use_validation.push_back(0);dology.push_back(0);add_ratio.push_back(1);rebin.push_back(1);

  hist.push_back("nElikeRing_angle_nring3_cut1_nring1_mulike0_michel0");
  input_type.push_back(3);mode_type.push_back(3);scale.push_back(1);
  use_validation.push_back(0);dology.push_back(0);add_ratio.push_back(1);rebin.push_back(1);

  hist.push_back("n_michel_electron_cut3_nring1_mulike0_michel0");
  input_type.push_back(3);mode_type.push_back(3);scale.push_back(1);
  use_validation.push_back(0);dology.push_back(0);add_ratio.push_back(1);rebin.push_back(1);

  hist.push_back("ntag_multiplicity_cut4_nring1_mulike0_michel0");
  input_type.push_back(3);mode_type.push_back(3);scale.push_back(1);
  use_validation.push_back(0);dology.push_back(0);add_ratio.push_back(1);rebin.push_back(1);*/

  //hist.push_back("mass_proton_reco_cut3_nring1_mulike0_michel0");
  //input_type.push_back(3);mode_type.push_back(3);scale.push_back(1);
  //use_validation.push_back(0);dology.push_back(1);add_ratio.push_back(1);rebin.push_back(5);

  hist.push_back("mom_proton_reco_cut3_nring1_mulike0_michel0");
  input_type.push_back(2);mode_type.push_back(2);scale.push_back(1);
  use_validation.push_back(0);dology.push_back(0);add_ratio.push_back(0);rebin.push_back(5);

  /*hist.push_back("true_angle_min_mid_lepton");
  input_type.push_back(3);mode_type.push_back(3);scale.push_back(1);
  use_validation.push_back(1);dology.push_back(0);add_ratio.push_back(1);rebin.push_back(1);

  hist.push_back("true_angle_min_max_lepton");
  input_type.push_back(3);mode_type.push_back(3);scale.push_back(1);
  use_validation.push_back(1);dology.push_back(0);add_ratio.push_back(1);rebin.push_back(1);

  hist.push_back("true_angle_mid_max_lepton");
  input_type.push_back(3);mode_type.push_back(3);scale.push_back(1);
  use_validation.push_back(1);dology.push_back(0);add_ratio.push_back(1);rebin.push_back(1);

  hist.push_back("angle_min_mom_lepton_proton_fp0");
  input_type.push_back(3);mode_type.push_back(3);scale.push_back(1);
  use_validation.push_back(1);dology.push_back(1);add_ratio.push_back(1);rebin.push_back(1);

  hist.push_back("angle_mid_mom_lepton_proton_fp0");
  input_type.push_back(3);mode_type.push_back(3);scale.push_back(1);
  use_validation.push_back(1);dology.push_back(1);add_ratio.push_back(1);rebin.push_back(1);

  hist.push_back("angle_max_mom_lepton_proton_fp0");
  input_type.push_back(3);mode_type.push_back(3);scale.push_back(1);
  use_validation.push_back(1);dology.push_back(1);add_ratio.push_back(1);rebin.push_back(1);

  hist.push_back("true_min_mom_lepton_fp0");
  input_type.push_back(3);mode_type.push_back(3);scale.push_back(1);
  use_validation.push_back(1);dology.push_back(1);add_ratio.push_back(1);rebin.push_back(1);

  hist.push_back("true_mid_mom_lepton_fp0");
  input_type.push_back(3);mode_type.push_back(3);scale.push_back(1);
  use_validation.push_back(1);dology.push_back(1);add_ratio.push_back(1);rebin.push_back(1);

  hist.push_back("true_max_mom_lepton_fp0");
  input_type.push_back(3);mode_type.push_back(3);scale.push_back(1);
  use_validation.push_back(1);dology.push_back(1);add_ratio.push_back(1);rebin.push_back(1);*/

  TFile *input_def,*input_up;
  TH1 *first_hist;
  for(int h=0;h<hist.size();h++){
    if(use_validation[h]) {
      input_def = TFile::Open(Form("../output/%s.%s.mode_%s_validation.root",type[input_type.at(h)].c_str(),period.c_str(),type[mode_type.at(h)].c_str()));
      input_up = TFile::Open(Form("../output/%s.%s.mode_%s_%s_validation.root",type[input_type.at(h)].c_str(),period.c_str(),type[mode_type.at(h)].c_str(),syst.c_str()));
    }else {
      input_def = TFile::Open(Form("../output/%s.%s.mode_%s.root",type[input_type.at(h)].c_str(),period.c_str(),type[mode_type.at(h)].c_str()));
      input_up = TFile::Open(Form("../output/%s.%s.mode_%s_%s.root",type[input_type.at(h)].c_str(),period.c_str(),type[mode_type.at(h)].c_str(),syst.c_str()));
    }
    TCanvas *c = new TCanvas(Form("canvas%d",h),"",800,600);
    if(add_ratio[h]){
      TPad* p1 = new TPad("main","main",0.0,0.2,1.0,1.0,10,0,0);
      p1->SetNumber(1);
      p1->SetBottomMargin(0);
      p1->Draw();
      TPad* p2 = new TPad("ratio","ratio",0.0,0.0,1.0,0.2,10,0,0);
      p2->SetTopMargin(0);
      p2->SetBottomMargin(0.5);
      p2->SetNumber(2);
      p2->Draw();
      c->cd(1);
      if(dology[h]) p1->SetLogy();
    }
    if(dology[h]) c->SetLogy();
    string save_name = "";
    save_name = "hist/compare_syst_" + syst;
    TH1* hist_def = (TH1*) input_def->Get(hist.at(h).c_str());
    TH1* hist_up = (TH1*) input_up->Get(hist.at(h).c_str());
    hist_def->Rebin(rebin[h]);
    hist_up->Rebin(rebin[h]);
    float default_integral_def = hist_def->Integral();
    if(scale[h]){
      hist_def->Scale(1./hist_def->Integral());
      hist_up->Scale(1./hist_up->Integral());
    }
    float maximum = hist_def->GetMaximum();
    if(hist_up->GetMaximum() > maximum) maximum = hist_up->GetMaximum();
    hist_def->SetLineWidth(2);
    hist_up->SetLineWidth(2);
    hist_up->SetLineColor(2);
    hist_def->GetYaxis()->SetRangeUser(0,maximum*1.2);
    if(dology[h]) {
      if(scale[h]) hist_up->GetYaxis()->SetRangeUser(1./default_integral_def, maximum*1.2);
      else hist_up->GetYaxis()->SetRangeUser(1, maximum*1.2);
    }
    hist_def->Draw("hist E0");
    hist_up->Draw("hist same E0");

    if(add_ratio[h]){
      c->cd(2);
      TH1 *ratio_up = (TH1*) hist_up->Clone("clone_hist_up");
      float xmin = ratio_up->GetBinLowEdge(1);
      float xmax = ratio_up->GetBinLowEdge(ratio_up->GetNbinsX())+ratio_up->GetBinWidth(ratio_up->GetNbinsX());
      //TH1* frame;
      //frame=gPad->DrawFrame(xmin, 0.5, xmax, 2.5);
      //frame->GetYaxis()->SetLabelSize(0.1);
      //frame->GetXaxis()->SetLabelSize(0.2);
      ratio_up->GetYaxis()->SetLabelSize(0.1);
      ratio_up->GetXaxis()->SetLabelSize(0.2);
      ratio_up->GetYaxis()->SetRangeUser(0.5,3);
      ratio_up->Divide(hist_def);
      ratio_up->Draw();
      TLine *line = new TLine(xmin,1,xmax,1);
      line->SetLineStyle(2);
      line->Draw();
      TLine *line2 = new TLine(xmin,2,xmax,2);
      line2->SetLineStyle(2);
      line2->Draw();
    }
    save_name += "_" + hist.at(h);
    save_name += "_" + type[input_type.at(h)] + "_" + type[mode_type.at(h)] ;
    save_name += "_" + period + ".pdf";
    c->SaveAs(save_name.c_str());

  }

}

#include "TGraphAsymmErrors.h"
void make_compare_plot_th1_syst(){
  string period = "sk4";
  string type[] = {"p_epi","p_mupi","p_eee_def","p_eee_miura","p_mumumu_def","p_mumumu_miura","p_emumu_def","p_emumu_miura","p_muee_def","p_muee_miura","p_eemu","fcmc","fcdt","single_mu_def","single_mu_miura"};
  string syst = "cd";

  gROOT->SetStyle("Plain");
  gStyle->SetOptStat(0);
  vector<string> hist;
  vector<int> dology,input_type,add_ratio,mode_type,use_validation,scale,rebin;
  hist.clear();dology.clear();input_type.clear();scale.clear();rebin.clear();
  add_ratio.clear();mode_type.clear();use_validation.clear();

  hist.push_back("total_true_mass_fp0");
  input_type.push_back(7);mode_type.push_back(7);scale.push_back(1);
  use_validation.push_back(1);dology.push_back(1);add_ratio.push_back(1);rebin.push_back(1);

  hist.push_back("total_true_mom_fp0");
  input_type.push_back(7);mode_type.push_back(7);scale.push_back(1);
  use_validation.push_back(1);dology.push_back(1);add_ratio.push_back(1);rebin.push_back(1);

  hist.push_back("true_mom_lepton_fp0");
  input_type.push_back(7);mode_type.push_back(7);scale.push_back(1);
  use_validation.push_back(1);dology.push_back(1);add_ratio.push_back(1);rebin.push_back(1);

  hist.push_back("nRing_cut1_nring1_mulike0_michel1");
  input_type.push_back(7);mode_type.push_back(7);scale.push_back(1);
  use_validation.push_back(0);dology.push_back(0);add_ratio.push_back(1);rebin.push_back(1);

  hist.push_back("nElikeRing_angle_nring3_cut1_nring1_mulike0_michel1");
  input_type.push_back(7);mode_type.push_back(7);scale.push_back(1);
  use_validation.push_back(0);dology.push_back(0);add_ratio.push_back(1);rebin.push_back(1);

  hist.push_back("n_michel_electron_cut3_nring1_mulike0_michel1");
  input_type.push_back(7);mode_type.push_back(7);scale.push_back(1);
  use_validation.push_back(0);dology.push_back(0);add_ratio.push_back(1);rebin.push_back(1);

  hist.push_back("ntag_multiplicity_cut4_nring1_mulike0_michel1");
  input_type.push_back(7);mode_type.push_back(7);scale.push_back(1);
  use_validation.push_back(0);dology.push_back(0);add_ratio.push_back(1);rebin.push_back(1);

  hist.push_back("mass_proton_reco_cut4_nring1_mulike0_michel1");
  input_type.push_back(7);mode_type.push_back(7);scale.push_back(1);
  use_validation.push_back(0);dology.push_back(1);add_ratio.push_back(1);rebin.push_back(5);

  hist.push_back("mom_proton_reco_cut4_nring1_mulike0_michel1");
  input_type.push_back(7);mode_type.push_back(7);scale.push_back(1);
  use_validation.push_back(0);dology.push_back(1);add_ratio.push_back(1);rebin.push_back(5);

  hist.push_back("true_angle_min_mid_lepton");
  input_type.push_back(7);mode_type.push_back(7);scale.push_back(1);
  use_validation.push_back(1);dology.push_back(0);add_ratio.push_back(1);rebin.push_back(1);

  hist.push_back("true_angle_min_max_lepton");
  input_type.push_back(7);mode_type.push_back(7);scale.push_back(1);
  use_validation.push_back(1);dology.push_back(0);add_ratio.push_back(1);rebin.push_back(1);

  hist.push_back("true_angle_mid_max_lepton");
  input_type.push_back(7);mode_type.push_back(7);scale.push_back(1);
  use_validation.push_back(1);dology.push_back(0);add_ratio.push_back(1);rebin.push_back(1);

  hist.push_back("angle_min_mom_lepton_proton_fp0");
  input_type.push_back(7);mode_type.push_back(7);scale.push_back(1);
  use_validation.push_back(1);dology.push_back(1);add_ratio.push_back(1);rebin.push_back(1);

  hist.push_back("angle_mid_mom_lepton_proton_fp0");
  input_type.push_back(7);mode_type.push_back(7);scale.push_back(1);
  use_validation.push_back(1);dology.push_back(1);add_ratio.push_back(1);rebin.push_back(1);

  hist.push_back("angle_max_mom_lepton_proton_fp0");
  input_type.push_back(7);mode_type.push_back(7);scale.push_back(1);
  use_validation.push_back(1);dology.push_back(1);add_ratio.push_back(1);rebin.push_back(1);

  hist.push_back("true_min_mom_lepton_fp0");
  input_type.push_back(7);mode_type.push_back(7);scale.push_back(1);
  use_validation.push_back(1);dology.push_back(1);add_ratio.push_back(1);rebin.push_back(1);

  hist.push_back("true_mid_mom_lepton_fp0");
  input_type.push_back(7);mode_type.push_back(7);scale.push_back(1);
  use_validation.push_back(1);dology.push_back(1);add_ratio.push_back(1);rebin.push_back(1);

  hist.push_back("true_max_mom_lepton_fp0");
  input_type.push_back(7);mode_type.push_back(7);scale.push_back(1);
  use_validation.push_back(1);dology.push_back(1);add_ratio.push_back(1);rebin.push_back(1);

  TFile *input_def,*input_down,*input_up;
  TH1 *first_hist;
  for(int h=0;h<hist.size();h++){
    if(use_validation[h]) {
      input_def = TFile::Open(Form("../output/%s.%s.mode_%s_validation.root",type[input_type.at(h)].c_str(),period.c_str(),type[mode_type.at(h)].c_str()));
      input_up = TFile::Open(Form("../output/%s.%s.mode_%s_%sup_validation.root",type[input_type.at(h)].c_str(),period.c_str(),type[mode_type.at(h)].c_str(),syst.c_str()));
      input_down = TFile::Open(Form("../output/%s.%s.mode_%s_%sdown_validation.root",type[input_type.at(h)].c_str(),period.c_str(),type[mode_type.at(h)].c_str(),syst.c_str()));
    }else {
      input_def = TFile::Open(Form("../output/%s.%s.mode_%s.root",type[input_type.at(h)].c_str(),period.c_str(),type[mode_type.at(h)].c_str()));
      input_up = TFile::Open(Form("../output/%s.%s.mode_%s_%sup.root",type[input_type.at(h)].c_str(),period.c_str(),type[mode_type.at(h)].c_str(),syst.c_str()));
      input_down = TFile::Open(Form("../output/%s.%s.mode_%s_%sdown.root",type[input_type.at(h)].c_str(),period.c_str(),type[mode_type.at(h)].c_str(),syst.c_str()));
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
    TH1* hist_down = (TH1*) input_down->Get(hist.at(h).c_str());
    hist_def->Rebin(rebin[h]);
    hist_up->Rebin(rebin[h]);
    hist_down->Rebin(rebin[h]);
    float default_integral_def = hist_def->Integral();
    if(scale[h]){
      hist_def->Scale(1./hist_def->Integral());
      hist_up->Scale(1./hist_up->Integral());
      hist_down->Scale(1./hist_down->Integral());
    }
    float maximum = hist_def->GetMaximum();
    if(hist_up->GetMaximum() > maximum) {
      maximum = hist_up->GetMaximum();
      if(hist_down->GetMaximum() > maximum) maximum = hist_down->GetMaximum();
    }
    else if (hist_down->GetMaximum() > maximum) maximum = hist_down->GetMaximum();
    hist_def->SetLineWidth(2);
    hist_up->SetLineWidth(2);
    hist_down->SetLineWidth(2);
    hist_up->SetLineColor(2);
    hist_down->SetLineColor(4);
    hist_up->GetYaxis()->SetRangeUser(0,maximum*1.2);
    if(dology[h]) {
      if(scale[h]) hist_up->GetYaxis()->SetRangeUser(1./default_integral_def, maximum*1.2);
      else hist_up->GetYaxis()->SetRangeUser(1, maximum*1.2);
    }
    hist_up->Draw("hist E0");
    hist_def->Draw("hist same E0");
    hist_down->Draw("hist same E0");

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
      TH1 *ratio_down = (TH1*) hist_down->Clone("clone_hist_down");
      ratio_down->GetYaxis()->SetRangeUser(0.5,2.5);
      ratio_down->Divide(hist_def);
      ratio_down->Draw("same");
      TLine *line = new TLine(xmin,1,xmax,1);
      line->SetLineStyle(2);
      line->Draw();
      TLine *line2 = new TLine(xmin,2,xmax,2);
      line2->SetLineStyle(2);
      line2->Draw();
    }
    /*if(scale[s] && this_hist->GetEntries()) this_hist->Scale(1./this_hist->Integral());
    this_hist->SetLineWidth(2);
    if(add_ratio[s]) c->cd(1);
    if(h==0){
      if(scale[s]) this_hist->GetYaxis()->SetRangeUser(0,evtmax_scale*1.2);
      else this_hist->GetYaxis()->SetRangeUser(0,evtmax*1.2);
      if(dology[s]){
        if(scale[s]) this_hist->GetYaxis()->SetRangeUser(evtmin_scale,evtmax_scale*1.2);
        else {
          this_hist->GetYaxis()->SetRangeUser(1,evtmax*1.2);
          if(very_small_evtmax<1) this_hist->GetYaxis()->SetRangeUser(very_small_evtmax*0.01,evtmax*1.2);
        }
      }
      this_hist->SetLineColor(1);
      this_hist->Draw("hist E0");
      //this_hist->Draw();
      first_hist = this_hist;
    }
    else{
      this_hist->SetLineColor(h+1);
      this_hist->Draw("same hist E0");
      TH1 *ratio_hist = (TH1*) this_hist->Clone("ratio_plot");
      ratio_hist->Divide(first_hist);
      float xmin = ratio_hist->GetBinLowEdge(1);
      float xmax = ratio_hist->GetBinLowEdge(ratio_hist->GetNbinsX())+ratio_hist->GetBinWidth(ratio_hist->GetNbinsX());
      if(add_ratio[s]){
        c->cd(2);
        if(h==1){
          TH1* frame;
          frame=gPad->DrawFrame(xmin, ratio_min*0.9, xmax, ratio_max*1.1);
          if(ratio_max<1) frame=gPad->DrawFrame(xmin, ratio_min*0.9, xmax, 1.1);
          if(ratio_min>1) frame=gPad->DrawFrame(xmin, 0.9, xmax, ratio_max*1.1);
          frame->GetYaxis()->SetLabelSize(0.1);
          frame->GetXaxis()->SetLabelSize(0.2);
        }
        ratio_hist->Draw("same");
        TLine *line = new TLine(xmin,1,xmax,1);
        line->SetLineStyle(2);
        line->Draw();
      }
    }*/
    save_name += "_" + hist.at(h);
    save_name += "_" + type[input_type.at(h)] + "_" + type[mode_type.at(h)] ;
    save_name += "_" + period + ".pdf";
    c->SaveAs(save_name.c_str());

  }

}

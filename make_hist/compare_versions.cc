#include <vector>
#include "TGraphAsymmErrors.h"
void compare_versions(){
  string type[] = {"p_epi","p_mupi","p_eee","p_mumumu","p_emumu","p_muee","p_eemu","fcmc","fcdt"};

  gROOT->SetStyle("Plain");
  gStyle->SetOptStat(0);
  vector<string> hist;
  vector<vector<string>> hist_set;
  vector<vector<int>> input_type_set,mode_type_set;
  vector<int> scale,dology,input_type,add_ratio,mode_type,use_validation;
  hist.clear();hist_set.clear();scale.clear();dology.clear();input_type.clear();
  input_type_set.clear();add_ratio.clear();mode_type_set.clear();mode_type.clear();use_validation.clear();

  hist.push_back("nRing_cut1_nring1_mulike0_michel1");
  hist.push_back("nRing_cut1_nring1_mulike0_michel1");
  input_type.push_back(7);input_type.push_back(7);
  mode_type.push_back(5);mode_type.push_back(5);
  add_ratio.push_back(1);scale.push_back(0);use_validation.push_back(0);
  hist_set.push_back(hist);dology.push_back(0);hist.clear();input_type_set.push_back(input_type);input_type.clear();
  mode_type_set.push_back(mode_type);mode_type.clear();

  hist.push_back("nMulikeRing_angle_nring3_cut1_nring1_mulike0_michel1");
  hist.push_back("nMulikeRing_angle_nring3_cut1_nring1_mulike0_michel1");
  input_type.push_back(7);input_type.push_back(7);
  mode_type.push_back(5);mode_type.push_back(5);
  add_ratio.push_back(1);scale.push_back(0);use_validation.push_back(0);
  hist_set.push_back(hist);dology.push_back(0);hist.clear();input_type_set.push_back(input_type);input_type.clear();
  mode_type_set.push_back(mode_type);mode_type.clear();

  hist.push_back("nMulikeRing_angle_nring2_cut1_nring1_mulike0_michel1");
  hist.push_back("nMulikeRing_angle_nring2_cut1_nring1_mulike0_michel1");
  input_type.push_back(7);input_type.push_back(7);
  mode_type.push_back(5);mode_type.push_back(5);
  add_ratio.push_back(1);scale.push_back(0);use_validation.push_back(0);
  hist_set.push_back(hist);dology.push_back(0);hist.clear();input_type_set.push_back(input_type);input_type.clear();
  mode_type_set.push_back(mode_type);mode_type.clear();

  hist.push_back("n_michel_electron_cut1_nring1_mulike0_michel1");
  hist.push_back("n_michel_electron_cut1_nring1_mulike0_michel1");
  input_type.push_back(7);input_type.push_back(7);
  mode_type.push_back(5);mode_type.push_back(5);
  add_ratio.push_back(1);scale.push_back(0);use_validation.push_back(0);
  hist_set.push_back(hist);dology.push_back(0);hist.clear();input_type_set.push_back(input_type);input_type.clear();
  mode_type_set.push_back(mode_type);mode_type.clear();

  hist.push_back("mass_proton_reco_cut3_nring1_mulike0_michel1");
  hist.push_back("mass_proton_reco_cut3_nring1_mulike0_michel1");
  input_type.push_back(7);input_type.push_back(7);
  mode_type.push_back(5);mode_type.push_back(5);
  add_ratio.push_back(1);scale.push_back(0);use_validation.push_back(0);
  hist_set.push_back(hist);dology.push_back(0);hist.clear();input_type_set.push_back(input_type);input_type.clear();
  mode_type_set.push_back(mode_type);mode_type.clear();

  hist.push_back("mom_proton_reco_cut3_nring1_mulike0_michel1");
  hist.push_back("mom_proton_reco_cut3_nring1_mulike0_michel1");
  input_type.push_back(7);input_type.push_back(7);
  mode_type.push_back(5);mode_type.push_back(5);
  add_ratio.push_back(1);scale.push_back(0);use_validation.push_back(0);
  hist_set.push_back(hist);dology.push_back(0);hist.clear();input_type_set.push_back(input_type);input_type.clear();
  mode_type_set.push_back(mode_type);mode_type.clear();

  TFile *input;
  TH1 *first_hist;
  for(int s=0;s<hist_set.size();s++){
    float evtmax=0,evtmax_scale=0,ratio_max=0,ratio_min=99999,evt_min,evtmin_scale=1;
    float very_small_evtmax = 99;
    for(int ss=0;ss<hist_set[s].size();ss++){//decide max event of the hist
      if(ss==0) input = TFile::Open(Form("../output/%s.sk4.mode_%s.root",type[input_type_set[s].at(ss)].c_str(),type[mode_type_set[s].at(ss)].c_str()));
      else input = TFile::Open(Form("../output/%s.sk4.mode_%s_new.root",type[input_type_set[s].at(ss)].c_str(),type[mode_type_set[s].at(ss)].c_str()));
      TH1* temp_hist = (TH1*) input->Get(hist_set[s].at(ss).c_str());
      if(temp_hist->GetMaximum() > evtmax) evtmax = temp_hist->GetMaximum();
      if(temp_hist->GetEntries() && temp_hist->GetMaximum()/temp_hist->Integral() > evtmax_scale){ 
        evtmax_scale = temp_hist->GetMaximum()/temp_hist->Integral();
        evtmin_scale = 1./temp_hist->Integral();
      }
      if(evtmax<1) very_small_evtmax = evtmax;
      if(ss==0) first_hist=temp_hist;
      else{
        TH1 *temp_ratio = (TH1*) temp_hist->Clone("temp_ratio");
        temp_ratio->Divide(first_hist);
        for(int n=0;n<temp_ratio->GetNbinsX();n++){
          float this_ratio = temp_ratio->GetBinContent(n+1);
          if(this_ratio<ratio_min) ratio_min=this_ratio;
          if(this_ratio>ratio_max) ratio_max=this_ratio;
        }
      }
    }
    TCanvas *c = new TCanvas(Form("canvas%d",s),"",800,600);
    if(add_ratio[s]){
      TPad* p1 = new TPad("main","main",0.0,0.2,1.0,1.0,10,0,0);
      p1->SetNumber(1);
      p1->SetBottomMargin(0);
      p1->Draw();
      TPad* p2 = new TPad("ratio","ratio",0.0,0.0,1.0,0.2,10,0,0);
      p2->SetTopMargin(0);
      p2->SetBottomMargin(0.5);
      p2->SetNumber(2);
      p2->Draw();
    }
    if(dology[s]) c->SetLogy();
    for(int h=0;h<hist_set[s].size();h++){
      if(h==0) input = TFile::Open(Form("../output/%s.sk4.mode_%s.root",type[input_type_set[s].at(h)].c_str(),type[mode_type_set[s].at(h)].c_str()));
      else input = TFile::Open(Form("../output/%s.sk4.mode_%s_new.root",type[input_type_set[s].at(h)].c_str(),type[mode_type_set[s].at(h)].c_str()));
      TH1* this_hist = (TH1*) input->Get(hist_set[s].at(h).c_str());
      if(scale[s] && this_hist->GetEntries()) this_hist->Scale(1./this_hist->Integral());
      this_hist->SetLineWidth(2);
      if(add_ratio[s]) c->cd(1);
      if(h==0){
        if(scale[s]) this_hist->GetYaxis()->SetRangeUser(0,evtmax_scale*1.2);
        else this_hist->GetYaxis()->SetRangeUser(0,evtmax*1.2);
        if(dology[s]){
          if(scale[s]) this_hist->GetYaxis()->SetRangeUser(evtmin_scale,evtmax_scale*1.2);
          else {
            this_hist->GetYaxis()->SetRangeUser(0.001,evtmax*1.2);
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
            //frame=gPad->DrawFrame(xmin, ratio_min*0.9, xmax, ratio_max*1.1);
            frame=gPad->DrawFrame(xmin, 0.7, xmax, 1.3);
            //if(ratio_max<1) frame=gPad->DrawFrame(xmin, ratio_min*0.9, xmax, 1.1);
            //if(ratio_min>1) frame=gPad->DrawFrame(xmin, 0.9, xmax, ratio_max*1.1);
            frame->GetYaxis()->SetLabelSize(0.1);
            frame->GetXaxis()->SetLabelSize(0.2);
          }
          ratio_hist->Draw("same");
          TLine *line = new TLine(xmin,1,xmax,1);
          line->SetLineStyle(2);
          line->Draw();
        }
      }
    }
    string save_name = "hist/compare_versions_";
    save_name += hist_set[s].at(0);
    save_name += "_" + type[input_type_set[s].at(0)] + "_" + type[mode_type_set[s].at(0)];
    save_name += ".pdf";
    cout << save_name << endl;
    c->SaveAs(save_name.c_str());
    //c->SaveAs("hist/temp.pdf");

  }

}

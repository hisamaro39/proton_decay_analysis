#include <vector>
#include "TGraphAsymmErrors.h"
void make_compare_plot_phd(){
  string period = "sk4";
  string type[] = {"p_epi","p_mupi","p_eee","p_eee_def","p_eee_miura",//4
    "p_eee_take","p_mumumu_def","p_mumumu_miura","p_mumumu_take","p_emumu_def",//9
    "p_emumu_miura","p_muee_def","p_muee_miura","p_eemu","fcmc",//14
    "fcdt","single_mu_def","single_mu_miura","subgev_multiring"};

  gROOT->SetStyle("Plain");
  gStyle->SetOptStat(0);
  vector<string> hist,save_name,hist_title,xtitle,ytitle,comment;
  vector<vector<string>> hist_set,hist_title_set,comment_set;
  vector<vector<int>> input_type_set,mode_type_set,color_set;
  vector<int> scale,dology,input_type,add_ratio,mode_type,use_validation,rebin,color,add_error;
  vector<float> xmax,xmin,ymax,ymin,xcomment,ycomment,legxmin,legymin,legxmax,legymax,step;
  hist.clear();hist_set.clear();scale.clear();dology.clear();input_type.clear();rebin.clear();save_name.clear();
  input_type_set.clear();add_ratio.clear();mode_type_set.clear();mode_type.clear();use_validation.clear();
  color_set.clear();color.clear();xmax.clear();xmin.clear();ymax.clear();ymin.clear();hist_title_set.clear();hist_title.clear();
  xtitle.clear();ytitle.clear();comment_set.clear();comment.clear();xcomment.clear();ycomment.clear();
  legxmin.clear();legymin.clear();legxmax.clear();legymax.clear();add_error.clear();step.clear();

  hist.push_back("true_min_mom_lepton_fp1");
  hist_title.push_back("Minimum");
  input_type.push_back(4);mode_type.push_back(4);color.push_back(2);
  hist.push_back("true_mid_mom_lepton_fp1");
  hist_title.push_back("Middle");
  input_type.push_back(4);mode_type.push_back(4);color.push_back(kGreen+2);
  hist.push_back("true_max_mom_lepton_fp1");
  hist_title.push_back("Maximum");
  input_type.push_back(4);mode_type.push_back(4);color.push_back(4);
  save_name.push_back("hist_phd/compare_true_mom_electrons_p_eee.pdf");
  xtitle.push_back("True momentum[MeV]");
  ytitle.push_back("Events");
  comment.push_back("SK4");
  comment.push_back("p#rightarrowe^{+}e^{+}e^{-}");
  comment.push_back("Free proton");
  xcomment.push_back(10);ycomment.push_back(240);
  step.push_back(15);
  ymin.push_back(0);ymax.push_back(250);xmin.push_back(0);xmax.push_back(500);
  legxmin.push_back(0.2);legymin.push_back(0.45);legxmax.push_back(0.5);legymax.push_back(0.75);
  add_ratio.push_back(0);
  add_error.push_back(0);
  scale.push_back(0);
  use_validation.push_back(1);
  rebin.push_back(1);
  dology.push_back(0);
  hist_set.push_back(hist);
  hist_title_set.push_back(hist_title);
  input_type_set.push_back(input_type);
  mode_type_set.push_back(mode_type);
  color_set.push_back(color);
  comment_set.push_back(comment);
  hist.clear();input_type.clear();mode_type.clear();color.clear();hist_title.clear();comment.clear();

  hist.push_back("total_true_mass_free");
  hist_title.push_back("Free proton");
  input_type.push_back(4);mode_type.push_back(4);color.push_back(1);
  hist.push_back("total_true_mass_sstate");
  hist_title.push_back("s-state");
  input_type.push_back(4);mode_type.push_back(4);color.push_back(2);
  hist.push_back("total_true_mass_pstate");
  hist_title.push_back("p-state");
  input_type.push_back(4);mode_type.push_back(4);color.push_back(kGreen+2);
  hist.push_back("total_true_mass_cor");
  hist_title.push_back("correlated decay");
  input_type.push_back(4);mode_type.push_back(4);color.push_back(4);
  save_name.push_back("hist_phd/compare_total_true_mass_hstate.pdf");
  xtitle.push_back("True proton mass[MeV]");
  ytitle.push_back("Events");
  comment.push_back("SK4");
  comment.push_back("p#rightarrowe^{+}e^{+}e^{-}");
  xcomment.push_back(10);ycomment.push_back(2000);
  step.push_back(800);
  ymin.push_back(1);ymax.push_back(3000);xmin.push_back(0);xmax.push_back(1000);
  legxmin.push_back(0.2);legymin.push_back(0.45);legxmax.push_back(0.5);legymax.push_back(0.75);
  add_ratio.push_back(0);
  add_error.push_back(0);
  scale.push_back(0);
  use_validation.push_back(1);
  rebin.push_back(1);
  dology.push_back(1);
  hist_set.push_back(hist);
  hist_title_set.push_back(hist_title);
  input_type_set.push_back(input_type);
  mode_type_set.push_back(mode_type);
  color_set.push_back(color);
  comment_set.push_back(comment);
  hist.clear();input_type.clear();mode_type.clear();color.clear();hist_title.clear();comment.clear();

  hist.push_back("total_true_mom_free");
  hist_title.push_back("Free proton");
  input_type.push_back(4);mode_type.push_back(4);color.push_back(1);
  hist.push_back("total_true_mom_sstate");
  hist_title.push_back("s-state");
  input_type.push_back(4);mode_type.push_back(4);color.push_back(2);
  hist.push_back("total_true_mom_pstate");
  hist_title.push_back("p-state");
  input_type.push_back(4);mode_type.push_back(4);color.push_back(kGreen+2);
  hist.push_back("total_true_mom_cor");
  hist_title.push_back("correlated decay");
  input_type.push_back(4);mode_type.push_back(4);color.push_back(4);
  save_name.push_back("hist_phd/compare_total_true_mom_hstate.pdf");
  xtitle.push_back("True proton momentum[MeV]");
  ytitle.push_back("Events");
  comment.push_back("SK4");
  comment.push_back("p#rightarrowe^{+}e^{+}e^{-}");
  xcomment.push_back(100);ycomment.push_back(2000);
  step.push_back(800);
  ymin.push_back(1);ymax.push_back(3000);xmin.push_back(0);xmax.push_back(800);
  legxmin.push_back(0.6);legymin.push_back(0.45);legxmax.push_back(0.9);legymax.push_back(0.75);
  add_ratio.push_back(0);
  add_error.push_back(0);
  scale.push_back(0);
  use_validation.push_back(1);
  rebin.push_back(1);
  dology.push_back(1);
  hist_set.push_back(hist);
  hist_title_set.push_back(hist_title);
  input_type_set.push_back(input_type);
  mode_type_set.push_back(mode_type);
  color_set.push_back(color);
  comment_set.push_back(comment);
  hist.clear();input_type.clear();mode_type.clear();color.clear();hist_title.clear();comment.clear();

  TFile *input;
  TH1 *first_hist;
  for(int s=0;s<hist_set.size();s++){
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
    TLegend *leg = new TLegend(legxmin[s],legymin[s],legxmax[s],legymax[s]);
    leg->SetFillColor(0);
    leg->SetBorderSize(0);
    for(int h=0;h<hist_set[s].size();h++){
      if(use_validation[s]) input = TFile::Open(Form("../output/%s.%s.mode_%s_validation.root",type[input_type_set[s].at(h)].c_str(),period.c_str(),type[mode_type_set[s].at(h)].c_str()));
      else input = TFile::Open(Form("../output/%s.%s.mode_%s.root",type[input_type_set[s].at(h)].c_str(),period.c_str(),type[mode_type_set[s].at(h)].c_str()));
      TH1* this_hist = (TH1*) input->Get(hist_set[s].at(h).c_str());
      this_hist->Rebin(rebin[s]);
      if(scale[s] && this_hist->GetEntries()) this_hist->Scale(1./this_hist->Integral());
      this_hist->SetLineWidth(2);
      this_hist->SetLineColor(color_set[s].at(h));
      leg->AddEntry(this_hist,hist_title_set[s].at(h).c_str(),"l");
      if(add_ratio[s]) c->cd(1);
      if(h==0){
        this_hist->GetYaxis()->SetRangeUser(ymin[s],ymax[s]);
        this_hist->GetXaxis()->SetRangeUser(xmin[s],xmax[s]);
        this_hist->GetXaxis()->SetTitle(xtitle[s].c_str());
        this_hist->GetYaxis()->SetTitle(ytitle[s].c_str());
        if(add_error[s]) this_hist->Draw("hist E0");
        else this_hist->Draw("hist");
        first_hist = this_hist;
      }
      else{
        if(add_error[s]) this_hist->Draw("same hist E0");
        else this_hist->Draw("same hist");
        TH1 *ratio_hist = (TH1*) this_hist->Clone("ratio_plot");
        ratio_hist->Divide(first_hist);
        if(add_ratio[s]){
          c->cd(2);
          if(h==1){
            TH1* frame;
            frame=gPad->DrawFrame(xmin[s], 0.5, xmax[s], 1.5);
            frame->GetYaxis()->SetLabelSize(0.1);
            frame->GetXaxis()->SetLabelSize(0.2);
          }
          ratio_hist->Draw("same");
          TLine *line = new TLine(xmin[s],1,xmax[s],1);
          line->SetLineStyle(2);
          line->Draw();
          TLine *line2 = new TLine(xmin[s],2,xmax[s],2);
          line2->SetLineStyle(2);
          line2->Draw();
        }
      }
    }
    leg->Draw();
    TLatex *text = new TLatex();
    text->SetTextSize(0.04);
    for(int t=0;t<comment_set[s].size();t++){
      cout << comment_set[s].at(t) << endl;
      text->DrawLatex(xcomment[s],ycomment[s]-step[s]*t,comment_set[s].at(t).c_str());
    }
    c->SaveAs(save_name[s].c_str());

  }

}

#include <vector>
void compare_data_mc(){
  string type[] = {"p_epi","p_mupi","p_eee","p_mumumu","p_emumu",//4
    "p_muee","fcmc","fcdt","subgev_multiring","subgev_onemulike"};
  int mode_id = 2;
  string sk_period = "sk4";
  bool add_signal = false;
  float signal_scale = 0.03;

  int range_momentum[] = {100,250,400,630,1000,2500,5000,10000,100000};
  int n_range_momentum = sizeof(range_momentum)/sizeof(int);

  gROOT->SetStyle("Plain");
  gStyle->SetOptStat(0);
  vector<string> hist;
  vector<vector<string>> hist_set;
  vector<vector<int>> input_type_set,mode_type_set;
  vector<int> scale,dology,input_type,mode_type,add_ratio,dorebin;
  hist.clear();hist_set.clear();scale.clear();dology.clear();dorebin.clear();input_type.clear();input_type_set.clear();add_ratio.clear();
  mode_type_set.clear();mode_type.clear();
  vector<string> hist_name;
  hist_name.clear();

  /*hist_name.push_back("distance_to_wall_thr50");
  add_ratio.push_back(1);scale.push_back(0);dology.push_back(0);dorebin.push_back(1);

  hist_name.push_back("distance_to_wall_scaled_thr50");
  add_ratio.push_back(1);scale.push_back(0);dology.push_back(0);dorebin.push_back(1);

  hist_name.push_back("distance_to_wall_thr100");
  add_ratio.push_back(1);scale.push_back(0);dology.push_back(0);dorebin.push_back(1);

  hist_name.push_back("distance_to_wall_scaled_thr100");
  add_ratio.push_back(1);scale.push_back(0);dology.push_back(0);dorebin.push_back(1);

  hist_name.push_back("distance_to_wall_thr150");
  add_ratio.push_back(1);scale.push_back(0);dology.push_back(0);dorebin.push_back(1);

  hist_name.push_back("distance_to_wall_scaled_thr150");
  add_ratio.push_back(1);scale.push_back(0);dology.push_back(0);dorebin.push_back(1);

  hist_name.push_back("n_michel_electron");
  add_ratio.push_back(1);scale.push_back(0);dology.push_back(0);dorebin.push_back(1);

  //hist_name.push_back("visible_energy_cut2_nring0_mulike0_michel0");
  //add_ratio.push_back(1);scale.push_back(0);dology.push_back(0);dorebin.push_back(2);

  hist_name.push_back("nhit_OD_cut1_nring0_mulike0_michel0");
  add_ratio.push_back(1);scale.push_back(0);dology.push_back(1);dorebin.push_back(2);
  
  hist_name.push_back("nRing_cut1_nring0_mulike0_michel0");
  add_ratio.push_back(1);scale.push_back(0);dology.push_back(1);dorebin.push_back(0);

  hist_name.push_back("nElikeRing_nring3_cut1_nring0_mulike0_michel0");
  add_ratio.push_back(1);scale.push_back(0);dology.push_back(1);dorebin.push_back(0);

  hist_name.push_back("nElikeRing_nring2_cut1_nring0_mulike0_michel0");
  add_ratio.push_back(1);scale.push_back(0);dology.push_back(1);dorebin.push_back(0);

  hist_name.push_back("nMulikeRing_nring3_cut1_nring0_mulike0_michel0");
  add_ratio.push_back(1);scale.push_back(0);dology.push_back(1);dorebin.push_back(0);

  hist_name.push_back("nMulikeRing_nring2_cut1_nring0_mulike0_michel0");
  add_ratio.push_back(1);scale.push_back(0);dology.push_back(1);dorebin.push_back(0);

  hist_name.push_back("n_michel_electron_cut3_nring0_mulike0_michel0");
  add_ratio.push_back(1);scale.push_back(0);dology.push_back(1);dorebin.push_back(0);

  hist_name.push_back("mass_pi0_reco_elike2_cut4_nring0_mulike0_michel0");
  add_ratio.push_back(1);scale.push_back(0);dology.push_back(1);dorebin.push_back(2);

  hist_name.push_back("ntag_multiplicity_cut5_nring0_mulike0_michel0");
  add_ratio.push_back(1);scale.push_back(0);dology.push_back(1);dorebin.push_back(0);*/

  hist_name.push_back("mom_proton_reco_cut6_nring0_mulike0_michel0");
  add_ratio.push_back(0);scale.push_back(0);dology.push_back(0);dorebin.push_back(5);

  hist_name.push_back("mass_proton_reco_cut6_nring0_mulike0_michel0");
  add_ratio.push_back(0);scale.push_back(0);dology.push_back(0);dorebin.push_back(5);

  TH1 *first_hist;
  TFile *input;
  for(int s=0;s<hist_name.size();s++){
    float evtmax=0,evtmax_scale=0,ratio_max=0,ratio_min=99999;
    for(int ss=0;ss<2;ss++){//decide max event of the hist
      if(ss==0) input = TFile::Open(Form("../output/%s.%s.mode_%s.root",type[7].c_str(),sk_period.c_str(),type[mode_id].c_str()));//data
      if(ss==1) input = TFile::Open(Form("../output/%s.%s.mode_%s.root",type[6].c_str(),sk_period.c_str(),type[mode_id].c_str()));//MC
      //TH1* temp_hist = (TH1*) input->Get(Form("%s_fp0",hist_name[s].c_str()));
      TH1* temp_hist = (TH1*) input->Get(hist_name[s].c_str());
      if(dorebin[s]) temp_hist->Rebin(dorebin[s]);
      if(temp_hist->GetMaximum() > evtmax) evtmax = temp_hist->GetMaximum();
      if(temp_hist->GetEntries() && temp_hist->GetMaximum()/temp_hist->GetEntries() > evtmax_scale) 
        evtmax_scale = temp_hist->GetMaximum()/temp_hist->GetEntries();
      if(ss==0) first_hist=temp_hist;
      else{
        TH1 *temp_ratio = (TH1*) temp_hist->Clone("temp_ratio");
        temp_ratio->Divide(first_hist);
        for(int n=0;n<temp_ratio->GetNbinsX();n++){
          float this_ratio = temp_ratio->GetBinContent(n+1);
          if(this_ratio>1e-5 && this_ratio<ratio_min) ratio_min=this_ratio;
          if(this_ratio>ratio_max) ratio_max=this_ratio;
        }
      }
    }
    
    //cout << "ratio min/max=" << ratio_min << "/" << ratio_max << endl;
    TCanvas *c = new TCanvas(Form("canvas%d",s),"",800,600);
    if(add_ratio[s]){
      TPad* p1 = new TPad("main","main",0.0,0.2,1.0,1.0,10,0,0);
      p1->SetNumber(1);
      p1->SetBottomMargin(0);
      if(dology[s]) p1->SetLogy();
      p1->Draw();
      TPad* p2 = new TPad("ratio","ratio",0.0,0.0,1.0,0.2,10,0,0);
      p2->SetTopMargin(0);
      p2->SetBottomMargin(0.5);
      p2->SetNumber(2);
      p2->Draw();
    }
    
    if(dology[s]) c->SetLogy();
    string save_name = "";
    save_name = "hist/compare_data_mc_";
    for(int h=0;h<4;h++){
      if(!add_signal && h==2) break;
      if(h==0) input = TFile::Open(Form("../output/%s.%s.mode_%s.root",type[7].c_str(),sk_period.c_str(),type[mode_id].c_str()));//data
      else if(h==1) input = TFile::Open(Form("../output/%s.%s.mode_%s.root",type[6].c_str(),sk_period.c_str(),type[mode_id].c_str()));//MC
      else input = TFile::Open(Form("../output/%s.%s.mode_%s.root",type[mode_id].c_str(),sk_period.c_str(),type[mode_id].c_str()));//signal
      //TH1 *this_hist = (TH1*) input->Get(Form("%s_fp0",hist_name[s].c_str()));
      TH1 *this_hist = (TH1*) input->Get(hist_name[s].c_str());
      if(h==3) this_hist = (TH1*) input->Get(Form("%s_fp1",hist_name[s].c_str()));
      if(dorebin[s]) this_hist->Rebin(dorebin[s]);
      if(scale[s] && this_hist->GetEntries()) this_hist->Scale(1./this_hist->GetEntries());
      this_hist->SetLineWidth(2);
      if(add_ratio[s]) c->cd(1);
      if(h==0){
        if(scale[s]) this_hist->GetYaxis()->SetRangeUser(0,evtmax_scale*1.2);
        else this_hist->GetYaxis()->SetRangeUser(0,evtmax*1.2);
        if(dology[s]) this_hist->GetYaxis()->SetRangeUser(1,evtmax*1.2);
        this_hist->SetLineStyle(0);
        this_hist->SetMarkerStyle(8);
        this_hist->Draw();
        first_hist = this_hist;
      }
      else{
        if(h==1) this_hist->SetLineColor(2);
        else {
          this_hist->Scale(signal_scale);
          this_hist->SetLineStyle(2);
          if(h==2) this_hist->SetLineColor(3);
          if(h==3) this_hist->SetLineColor(6);
        }
        this_hist->Draw("same hist");
        if(h>1) continue;
        TH1 *ratio_hist = (TH1*) first_hist->Clone("ratio_plot");
        ratio_hist->Divide(this_hist);
        float xmin = ratio_hist->GetBinLowEdge(1);
        float xmax = ratio_hist->GetBinLowEdge(ratio_hist->GetNbinsX())+ratio_hist->GetBinWidth(ratio_hist->GetNbinsX());
        if(add_ratio[s]){
          c->cd(2);
          TH1* frame;
          frame=gPad->DrawFrame(xmin, ratio_min*0.9, xmax, ratio_max*1.1);
          if(ratio_max<1) frame=gPad->DrawFrame(xmin, ratio_min*0.9, xmax, 1.1);
          if(ratio_min>1) frame=gPad->DrawFrame(xmin, 0.9, xmax, ratio_max*1.1);
          frame->GetYaxis()->SetLabelSize(0.1);
          frame->GetXaxis()->SetLabelSize(0.2);
          TLine *line = new TLine(xmin,1,xmax,1);
          line->SetLineStyle(2);
          line->SetLineWidth(2);
          line->SetLineColor(1);
          line->Draw();
          ratio_hist->Draw("same");
        }
      }
    }
    save_name += hist_name[s] +  "_mode_" + type[mode_id] + "_"; 
    if(dology[s]) save_name += "logy_"; 
    if(add_signal) save_name += "add_signal_";
    save_name += sk_period + ".pdf";
    cout << "save_name = " << save_name << endl;
    c->SaveAs(save_name.c_str());
  }

}

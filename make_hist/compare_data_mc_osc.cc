#include <vector>
void compare_data_mc_osc(){
  string type[] = {"p_epi","p_mupi","p_eee","p_mumumu","p_emumu","p_muee","fcmc","fcdt"};
  bool add_signal = true;

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

  for(int t=1;t<=14;t++){
    for(int r=0;r<n_range_momentum-1;r++){
      hist.push_back(Form("zenith_angle_type%d_mom%d_%d",t,range_momentum[r],range_momentum[r+1]));
      hist.push_back(Form("zenith_angle_type%d_mom%d_%d",t,range_momentum[r],range_momentum[r+1]));
      hist.push_back(Form("zenith_angle_type%d_mom%d_%d_osc",t,range_momentum[r],range_momentum[r+1]));
      hist_set.push_back(hist);hist.clear();
      input_type.push_back(7);input_type.push_back(6);input_type.push_back(6);
      input_type_set.push_back(input_type);input_type.clear();
      mode_type.push_back(6);mode_type.push_back(6);mode_type.push_back(6);
      mode_type_set.push_back(mode_type);mode_type.clear();
      add_ratio.push_back(1);scale.push_back(0);dology.push_back(0);dorebin.push_back(0);
    }
  }

  TH1 *first_hist;
  TFile *input;
  for(int s=0;s<hist_set.size();s++){
    float evtmax=0,evtmax_scale=0,ratio_max=0,ratio_min=99999;
    for(int ss=0;ss<hist_set[s].size();ss++){//decide max event of the hist
      if(ss==3) continue;
      input = TFile::Open(Form("../output/%s.sk4.mode_%s.root",type[input_type_set[s].at(ss)].c_str(),type[mode_type_set[s].at(ss)].c_str()));
      TH1* temp_hist = (TH1*) input->Get(hist_set[s].at(ss).c_str());
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
    save_name = "hist_osc/compare_data_mc_";
    for(int h=0;h<hist_set[s].size();h++){
      input = TFile::Open(Form("../output/%s.sk4.mode_%s.root",type[input_type_set[s].at(h)].c_str(),type[mode_type_set[s].at(h)].c_str()));
      TH1 *this_hist = (TH1*) input->Get(hist_set[s].at(h).c_str());
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
        if(h==2) this_hist->SetLineColor(4);
        this_hist->Draw("same hist");
        if(h==3) {
          this_hist->Scale(0.01);
          this_hist->SetLineColor(6);
          this_hist->SetLineStyle(2);
          continue;
        }
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
            TLine *line = new TLine(xmin,1,xmax,1);
            line->SetLineStyle(2);
            line->SetLineWidth(2);
            line->SetLineColor(1);
          line->Draw();
          }
          ratio_hist->Draw("same");
        }
      }
      save_name += hist_set[s].at(h);
      save_name += "_input_" + type[input_type_set[s].at(h)]; 
      save_name += "_mode_" + type[mode_type_set[s].at(h)] + "_"; 
    }
    if(dology[s]) save_name += "logy_"; 
    if(dorebin[s]) {
      stringstream str;
      str << dorebin[s];
      save_name += "rebin" + str.str() + "_";
    }
    if(add_signal) save_name += "add_signal";
    save_name += ".pdf";
    c->SaveAs(save_name.c_str());

  }

}

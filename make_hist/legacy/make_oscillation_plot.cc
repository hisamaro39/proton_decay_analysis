#include <vector>
void make_oscillation_plot(){
  TFile *input_mc = TFile::Open("../output/fcmc.sk4.mode_p_epi.root");
  TFile *input_data = TFile::Open("../output/fcdt.sk4.mode_p_epi.root");
  string input_name[] = {"p_epi","p_mupi","p_eee","p_mumumu","p_emumu","p_muee","fcmc","fcdt"};

  int range_momentum[] = {100,250,400,630,1000,2500,5000,10000,100000};
  int n_range_momentum = sizeof(range_momentum)/sizeof(int);

  gROOT->SetStyle("Plain");
  gStyle->SetOptStat(0);
  vector<string> hist;
  vector<vector<string>> hist_set;
  vector<vector<int>> input_type_set;
  vector<int> scale,dology,input_type,add_ratio;
  hist.clear();hist_set.clear();scale.clear();dology.clear();input_type.clear();input_type_set.clear();add_ratio.clear();

  TH1 *data_hist,*mc_hist[2];
  for(int t=1;t<2;t++){
    //for(int r=0;r<n_range_momentum-1;r++){
    for(int r=0;r<1;r++){
      data_hist = (TH1*) input_data->Get(Form("zenith_angle_type%d_mom%d_%d",t,range_momentum[r],range_momentum[r+1]));
      mc_hist[0] = (TH1*) input_mc->Get(Form("zenith_angle_type%d_mom%d_%d",t,range_momentum[r],range_momentum[r+1]));
      mc_hist[1] = (TH1*) input_mc->Get(Form("zenith_angle_type%d_mom%d_%d_osc",t,range_momentum[r],range_momentum[r+1]));
      data_hist->GetYaxis()->SetRangeUser(0,300);
      data_hist->Draw();
      mc_hist[0]->SetLineColor(2);
      mc_hist[0]->Draw("hist same");
      mc_hist[1]->SetLineColor(4);
      mc_hist[1]->Draw("hist same");
    }
  }

  /*
  //hist list
  hist.push_back("zenith_angle_type1_mom100_250");hist.push_back("zenith_angle_type1_mom100_250");
  hist.push_back("zenith_angle_type1_mom100_250_osc");
  add_ratio.push_back(1);scale.push_back(0);input_type.push_back(7);input_type.push_back(6);input_type.push_back(6);
  hist_set.push_back(hist);dology.push_back(0);hist.clear();input_type_set.push_back(input_type);input_type.clear();

  TH1 *first_hist;
  for(int s=0;s<hist_set.size();s++){
    float evtmax=0,evtmax_scale=0,ratio_max=0,ratio_min=99999;
    for(int ss=0;ss<hist_set[s].size();ss++){//decide max event of the hist
      TH1* temp_hist = (TH1*) input[input_type_set[s].at(ss)]->Get(hist_set[s].at(ss).c_str());
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
      p1->Draw();
      TPad* p2 = new TPad("ratio","ratio",0.0,0.0,1.0,0.2,10,0,0);
      p2->SetTopMargin(0);
      p2->SetBottomMargin(0.5);
      p2->SetNumber(2);
      p2->Draw();
    }
    if(dology[s]) c->SetLogy();
    string save_name = "";
    save_name = "hist/compare_";
    for(int h=0;h<hist_set[s].size();h++){
      TH1* this_hist = (TH1*) input[input_type_set[s].at(h)]->Get(hist_set[s].at(h).c_str());
      if(scale[s] && this_hist->GetEntries()) this_hist->Scale(1./this_hist->GetEntries());
      this_hist->SetLineWidth(2);
      if(add_ratio[s]) c->cd(1);
      if(h==0){
        if(scale[s]) this_hist->GetYaxis()->SetRangeUser(0,evtmax_scale*1.2);
        else this_hist->GetYaxis()->SetRangeUser(0,evtmax*1.2);
        if(dology[s]) this_hist->GetYaxis()->SetRangeUser(1,evtmax*1.2);
        this_hist->SetLineColor(1);
        this_hist->Draw("hist");
        first_hist = this_hist;
      }
      else{
        this_hist->SetLineColor(h+1);
        this_hist->Draw("same hist");
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
      }
      save_name += hist_set[s].at(h);
      save_name += "_" + input_name[input_type_set[s].at(h)] + "_";
    }
    save_name += ".pdf";
    c->SaveAs(save_name.c_str());

  }
*/

}

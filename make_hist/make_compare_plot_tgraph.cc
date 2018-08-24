#include <vector>
void make_compare_plot_tgraph(){
  TFile *input[2];
  input[0] = TFile::Open("output/efficiency.root");

  gROOT->SetStyle("Plain");
  gStyle->SetOptStat(0);
  vector<string> hist;
  vector<vector<string>> hist_set;
  vector<vector<int>> input_type_set;
  vector<int> input_type;
  hist.clear();hist_set.clear();input_type.clear();input_type_set.clear();
  //hist list
  hist.push_back("efficiency_mom_gamma_vector_pid_elike_gamma_match_ring_only_epi0");
  hist.push_back("efficiency_mom_gamma_vector_prob_ring_angle_elike_gamma_match_ring_only_epi0");
  hist.push_back("efficiency_mom_gamma_vector_prob_hit_elike_gamma_match_ring_only_epi0");
  input_type.push_back(0);input_type.push_back(0);input_type.push_back(0);
  hist_set.push_back(hist);hist.clear();input_type_set.push_back(input_type);input_type.clear();
  hist.push_back("efficiency_mom_e_vector_pid_elike_e_match_ring_only_eee");
  hist.push_back("efficiency_mom_e_vector_prob_ring_angle_elike_e_match_ring_only_eee");
  hist.push_back("efficiency_mom_e_vector_prob_hit_elike_e_match_ring_only_eee");
  input_type.push_back(0);input_type.push_back(0);input_type.push_back(0);
  hist_set.push_back(hist);hist.clear();input_type_set.push_back(input_type);input_type.clear();
  hist.push_back("efficiency_mom_mu_vector_pid_mulike_mu_match_ring_only_mumumu");
  hist.push_back("efficiency_mom_mu_vector_prob_ring_angle_mulike_mu_match_ring_only_mumumu");
  hist.push_back("efficiency_mom_mu_vector_prob_hit_mulike_mu_match_ring_only_mumumu");
  input_type.push_back(0);input_type.push_back(0);input_type.push_back(0);
  hist_set.push_back(hist);hist.clear();input_type_set.push_back(input_type);input_type.clear();

  for(int s=0;s<hist_set.size();s++){
  //for(int s=0;s<1;s++){
    TCanvas *c = new TCanvas(Form("canvas%d",s),"",800,600);
    TPad* p1 = new TPad("main","main",0.0,0.2,1.0,1.0,10,0,0);
    p1->SetNumber(1);
    p1->SetBottomMargin(0);
    p1->Draw();
    TPad* p2 = new TPad("ratio","ratio",0.0,0.0,1.0,0.2,10,0,0);
    p2->SetTopMargin(0);
    p2->SetBottomMargin(0.5);
    p2->SetNumber(2);
    p2->Draw();
    string save_name = "";
    save_name = "hist/compare_";
    TGraphAsymmErrors *first_graph;
    double xmin=0,xmax=0;
    for(int h=0;h<hist_set[s].size();h++){
      TGraphErrors *ratio_graph = new TGraphErrors();
      TGraphAsymmErrors* this_graph = (TGraphAsymmErrors*) input[input_type_set[s].at(h)]->Get(hist_set[s].at(h).c_str());
      this_graph->SetLineWidth(2);
      if(h==0){
        int nbins = this_graph->GetN();
        double *x = this_graph->GetX();
        double *exhigh = this_graph->GetEXhigh();
        double *exlow = this_graph->GetEXlow();
        xmin = x[0] - exlow[0];
        xmax = x[nbins-1] + exhigh[nbins-1];
        c->cd(1);
        TH1* frame=gPad->DrawFrame(xmin, 0, xmax, 1.1);
        this_graph->SetLineColor(1);
        this_graph->GetYaxis()->SetRangeUser(0,1.1);
        this_graph->Draw("p");
        first_graph = this_graph;
      }
      else{
        c->cd(1);
        this_graph->SetLineColor(h+1);
        this_graph->Draw("p");
        double *x = first_graph->GetX();
        double *exlow = first_graph->GetEXlow();
        double *y1 = first_graph->GetY();
        double *eyhigh1 = first_graph->GetEYhigh();
        double *eylow1 = first_graph->GetEYlow();
        double *y2 = this_graph->GetY();
        double *eyhigh2 = this_graph->GetEYhigh();
        double *eylow2 = this_graph->GetEYlow();
        for(int b=0;b<first_graph->GetN();b++){
          float ratio = (y1[b]>0)? y2[b]/y1[b] : 0;
          float err = (y1[b]>0)? sqrt(y1[b]*y1[b]*eylow2[b]*eylow2[b]+y2[b]*y2[b]*eylow1[b]*eylow1[b])/pow(y1[b],2) : 0;
          //cout << "bin" << b << endl;
          //cout << "ratio=" << ratio << " +- " << err << endl;
          ratio_graph->SetPoint(b+1,x[b],ratio);
          ratio_graph->SetPointError(b+1,exlow[b],err);
        }
        c->cd(2);
        if(h==1) {
          TH1* frame=gPad->DrawFrame(xmin, 0.8, xmax, 1.2);
          frame->GetYaxis()->SetLabelSize(0.1);
          frame->GetXaxis()->SetLabelSize(0.2);
        }
        ratio_graph->SetLineWidth(2);
        ratio_graph->SetLineColor(h+1);
        ratio_graph->Draw("p");
      }
      save_name += hist_set[s].at(h);
      save_name += "_";
    }
    save_name += ".pdf";
    c->SaveAs(save_name.c_str());

  }

}

#include <vector>
void make_single_plot_tgraph(){
  TFile *input[2];
  input[0] = TFile::Open("output/efficiency.root");

  gROOT->SetStyle("Plain");
  gStyle->SetOptStat(0);
  vector<string> hist;
  vector<int> input_type;
  hist.clear();input_type.clear();
  //hist list
  //hist.push_back("efficiency_true_mom_lepton_match_ring_nring3_p_eee");
  //input_type.push_back(0);
  //hist.push_back("efficiency_true_mom_lepton_match_ring_nring3_p_mumumu");
  //input_type.push_back(0);
  //hist.push_back("efficiency_true_mom_muon_match_ring_nring3_p_muee");
  //input_type.push_back(0);
  //hist.push_back("efficiency_true_mom_lepton_match_ring_angle_elike_nring3_p_eee");
  //input_type.push_back(0);
  hist.push_back("efficiency_true_mom_muon_match_ring_angle_mulike_nring3_p_muee");
  input_type.push_back(0);

  for(int h=0;h<hist.size();h++){
    TCanvas *c = new TCanvas(Form("canvas%d",h),"",800,600);
    TGraphAsymmErrors* this_graph = (TGraphAsymmErrors*) input[input_type[h]]->Get(hist[h].c_str());
    this_graph->SetLineWidth(2);
    int nbins = this_graph->GetN();
    double *x = this_graph->GetX();
    double *exhigh = this_graph->GetEXhigh();
    double *exlow = this_graph->GetEXlow();
    xmin = x[0] - exlow[0];
    xmax = x[nbins-1] + exhigh[nbins-1];
    TH1* frame=gPad->DrawFrame(xmin, 0, xmax, 1.1);
    this_graph->SetLineColor(1);
    this_graph->GetYaxis()->SetRangeUser(0,1.1);
    this_graph->Draw("p");
    c->SaveAs(Form("hist/single_%s.pdf",hist[h].c_str()));
  }

}

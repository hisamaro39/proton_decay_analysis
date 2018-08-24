#include <vector>
void make_single_plot_th1(){
  TFile *input[6];
  input[0] = TFile::Open("../output/p_epi.sk4.mode_p_epi.root");
  input[1] = TFile::Open("../output/p_mupi.sk4.mode_p_mupi.root");
  input[2] = TFile::Open("../output/p_eee.sk4.mode_p_eee.root");
  input[3] = TFile::Open("../output/p_mumumu.sk4.mode_p_mumumu.root");
  input[4] = TFile::Open("../output/p_emumu.sk4.mode_p_emumu.root");
  input[5] = TFile::Open("../output/p_muee.sk4.mode_p_muee.root");
  string input_name[] = {"p_epi","p_mupi","p_eee","p_mumumu","p_emumu","p_muee"};

  gROOT->SetStyle("Plain");
  gStyle->SetOptStat(0);
  vector<string> hist;
  vector<int> dology,input_type,make_cutflow;
  hist.clear();dology.clear();input_type.clear();make_cutflow.clear();
  //hist list
  //hist.push_back("nRing_cut1_nring0_mulike0_michel0");  
  //input_type.push_back(2);dology.push_back(0);
  //hist.push_back("nElikeRing_nring3_cut1_nring0_mulike0_michel0");  
  //input_type.push_back(2);dology.push_back(0);
  //hist.push_back("nElikeRing_nring2_cut1_nring0_mulike0_michel0");  
  //input_type.push_back(2);dology.push_back(0);
  hist.push_back("cut_flow_nring0_mulike0_michel0");  
  input_type.push_back(0);dology.push_back(0);make_cutflow.push_back(1);

  for(int h=0;h<hist.size();h++){
    TCanvas *c = new TCanvas(Form("canvas%d",h),"",800,600);
    if(dology[h]) c->SetLogy();
    string save_name = "";
    TH1* this_hist = (TH1*) input[input_type[h]]->Get(hist[h].c_str());
    //this_hist->Rebin(2);
    this_hist->SetLineWidth(2);
    this_hist->SetLineColor(1);
    if(make_cutflow[h]) this_hist->Scale(1./this_hist->GetBinContent(1));
    this_hist->Draw("hist");
    c->SaveAs(Form("hist/single_%s_%s.pdf",hist[h].c_str(),input_name[input_type[h]].c_str()));
  }

}

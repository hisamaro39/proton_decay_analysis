#include <vector>
void open_signal_box(){
  string input_type = "fcdt_final";
  //string input_type = "fcmc_real";
  string mode_type = "p_muee";
  int period = 1;
  int nring = 1;
  int mulike = 0;
  int michel = 1;

  gROOT->SetStyle("Plain");
  gStyle->SetOptStat(0);
  gStyle->SetPalette(1);

  TFile *input;
  if(period==5) input = TFile::Open(Form("../output/%s.sk1_4.mode_%s.root",input_type.c_str(),mode_type.c_str()));
  else input = TFile::Open(Form("../output/%s.sk%d.mode_%s.root",input_type.c_str(),period,mode_type.c_str()));
  TCanvas *c = new TCanvas("canvas","",800,600);
  TH1* frame=gPad->DrawFrame(0, 0, 1250, 1000);
  string save_name = "";
  TGraph* graph = (TGraph*) input->Get(Form("mass_mom_proton_reco_ntag_cut5_nring%d_mulike%d_michel%d",nring,mulike,michel));  
  graph->SetMarkerStyle(7);
  TBox *box = new TBox(800,0,1050,100);
  box->SetFillStyle(0);
  box->SetLineWidth(2);
  box->SetLineColor(1);
  box->Draw();
  TBox *box2 = new TBox(800,100,1050,250);
  box2->SetFillStyle(0);
  box2->SetLineWidth(2);
  box2->SetLineColor(1);
  box2->Draw();
  graph->Draw("p");
  //c->SaveAs(Form("hist/single_mass_mom_proton_input_%s_mode_%s.pdf",type[input_type].c_str(),type[mode_type].c_str()));

}

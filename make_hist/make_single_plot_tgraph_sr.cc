#include <vector>
void make_single_plot_tgraph_sr(){
  string type[] = {"p_epi","p_mupi","p_eee","p_muee","p_mumumu","p_emumu","p_mumue","p_muee","fcmc","fcdt"};
  int input_type = 8;
  int mode_type = 2;
  int nring = 1;
  int mulike = 0;
  int michel = 0;

  gROOT->SetStyle("Plain");
  gStyle->SetOptStat(0);
  gStyle->SetPalette(1);
  //vector<string> graph;
  //vector<int> input_type,mode_type;
  //graph.clear();input_type.clear();mode_type.clear();
  //hist list
  //input_type.push_back(8);mode_type.push_back(2);

  TFile *input;
  input = TFile::Open(Form("../output/%s.sk4.mode_%s.root",type[input_type].c_str(),type[mode_type].c_str()));
  TCanvas *c = new TCanvas("canvas","",800,600);
  TH1* frame=gPad->DrawFrame(0, 0, 1250, 1000);
  string save_name = "";
  TGraph* graph_all = (TGraph*) input->Get(Form("mass_mom_proton_reco_cut5_nring%d_mulike%d_michel%d",nring,mulike,michel));  
  TGraph* graph_sr_low = (TGraph*) input->Get(Form("mass_mom_proton_reco_cut6_nring%d_mulike%d_michel%d",nring,mulike,michel));  
  TGraph* graph_sr_high = (TGraph*) input->Get(Form("mass_mom_proton_reco_cut7_nring%d_mulike%d_michel%d",nring,mulike,michel));  
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
  TEllipse *e = new TEllipse(925, 125, 125, 125);
  e->SetFillStyle(0);
  e->SetLineColor(2);
  e->SetLineWidth(2);
  e->Draw();
  TEllipse *e2 = new TEllipse(925, 125, 500, 500);
  e2->SetFillStyle(0);
  e2->SetLineColor(2);
  e2->SetLineWidth(2);
  e2->SetLineStyle(2);
  e2->Draw();
  TEllipse *e3 = new TEllipse(925, 125, 1000, 1000);
  e3->SetFillStyle(0);
  e3->SetLineColor(2);
  e3->SetLineWidth(2);
  e3->SetLineStyle(2);
  e3->Draw();
  graph_all->Draw("p");
  graph_sr_low->SetMarkerStyle(8);
  graph_sr_low->Draw("p");
  graph_sr_high->SetMarkerStyle(8);
  graph_sr_high->Draw("p");
  c->SaveAs(Form("hist/single_mass_mom_proton_input_%s_mode_%s.pdf",type[input_type].c_str(),type[mode_type].c_str()));

}

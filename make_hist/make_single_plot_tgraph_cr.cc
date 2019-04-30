#include <vector>
void make_single_plot_tgraph_cr(){
  string input_type = "fcdt_final";
  string mode_type = "p_mumumu";

  gROOT->SetStyle("Plain");
  gStyle->SetOptStat(0);
  gStyle->SetPalette(1);
  vector<string> graph;
  graph.clear();
  //graph list
  graph.push_back("mass_mom_proton_reco_cut5_nring0_mulike0_michel3");  

  TFile *input;
  TBox *box = new TBox(800,0,1050,100);
  box->SetFillStyle(0);
  box->SetLineWidth(2);
  TBox *box2 = new TBox(800,100,1050,250);
  box2->SetFillStyle(0);
  box2->SetLineWidth(2);
  TBox *box3 = new TBox(600,0,1250,450);
  box3->SetFillStyle(0);
  box3->SetLineWidth(2);
  box3->SetLineColor(2);
  for(int p=1;p<6;p++){
    if(p==5) input = TFile::Open(Form("../output/%s.sk1_4.mode_%s_CR.root",input_type.c_str(),mode_type.c_str()));
    else input = TFile::Open(Form("../output/%s.sk%d.mode_%s_CR.root",input_type.c_str(),p,mode_type.c_str()));
    //cout << Form("../output/%s.sk%d.mode_%s_CR.root",input_type.c_str(),p,mode_type.c_str()) << endl;
    //for(int h=0;h<graph.size();h++){
      TCanvas *c = new TCanvas(Form("canvas%d",p),"",800,600);
      TH1* frame=gPad->DrawFrame(0, 0, 2000, 2000);
      box->Draw();
      box2->Draw();
      box3->Draw();
      //string save_name = "";
      TGraph* this_graph = (TGraph*) input->Get(graph[0].c_str());
      this_graph->SetMarkerStyle(7);
      this_graph->Draw("p");
      //c->SaveAs(Form("hist/single_CR_%s_input_%s_mode_%s_sk%d.pdf",graph[h].c_str(),input_type.c_str(),mode_type.c_str(),p));
    //}
  }

}

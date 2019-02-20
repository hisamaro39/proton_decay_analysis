void test2(){

  TFile *input = TFile::Open("../output/fcmc.sk4.mode_p_eee.root");
  gROOT->SetStyle("Plain");
  gStyle->SetOptStat(0);
  TCanvas *c = new TCanvas("c","",800,600);
  TH1* frame=gPad->DrawFrame(0, 0, 1250, 1000);
  TGraph* graph = (TGraph*) input->Get("mass_mom_proton_reco_cut5_nring1_mulike0_michel0");
  graph->Draw("p");
  TGraph* graph2 = (TGraph*) input->Get("mass_mom_proton_reco_cut6_nring1_mulike0_michel0");
  graph2->SetMarkerStyle(8);
  graph2->Draw("p");
  TGraph* graph3 = (TGraph*) input->Get("mass_mom_proton_reco_cut7_nring1_mulike0_michel0");
  graph3->SetMarkerStyle(8);
  graph3->Draw("p");

  TBox *box = new TBox(800,0,1050,100);
  box->SetFillStyle(0);
  box->SetLineWidth(2);
  box->SetLineColor(1);
  TBox *box2 = new TBox(800,100,1050,250);
  box2->SetFillStyle(0);
  box2->SetLineWidth(2);
  box2->SetLineColor(1);

  box->Draw();
  box2->Draw();

  c->SaveAs("hist/mass_mom_eee_bkg.pdf");

}

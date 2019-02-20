void cherenkov(){

  gROOT->SetStyle("Plain");
  gStyle->SetOptStat(0);
  gStyle->SetPadBottomMargin(0.1);

  float n = 1.333;
  float mass_e = 0.511;//MeV
  float mass_mu = 105.7;//MeV
  TGraph *g_mu = new TGraph();
  TGraph *g_e = new TGraph();
  for(int p=120;p<10000;p++){
    float mom = p;//MeV
    float theta_mu = acos(sqrt(1+pow(mass_mu,2)/pow(mom,2))/n);
    float angle_mu = 180*theta_mu/3.14159;
    //cout << "muon mom[MeV]/theta/angle=" << mom << "/" << theta_mu << "/" << angle_mu << endl;
    g_mu->SetPoint(p-120,mom,angle_mu);
  }
  for(int p=60;p<1000000;p++){
    float mom = 0.01*p;//MeV
    float theta_e = acos(sqrt(1+pow(mass_e,2)/pow(mom,2))/n);
    float angle_e = 180*theta_e/3.14159;
    //cout << "electron mom[MeV]/thetai/angle=" << mom << "/" << theta_e << "/" << angle_e << endl;
    g_e->SetPoint(p-60,mom,angle_e);
  }
  TCanvas *c1 = new TCanvas("c1","",800,600);
  TH1* frame = gPad->DrawFrame(0.1, 0, 10000, 50);
  frame->GetXaxis()->SetLabelSize(0.04);
  frame->GetXaxis()->SetTitle("Momentum [MeV]");
  frame->GetXaxis()->SetTitleOffset(1.1);
  frame->GetYaxis()->SetLabelSize(0.04);
  frame->GetYaxis()->SetTitle("Opening angle [degree]");
  c1->SetLogx();
  g_mu->SetLineWidth(2);
  g_mu->SetLineColor(2);
  g_e->SetLineWidth(2);
  g_e->SetLineColor(4);
  g_mu->Draw("l");
  g_e->Draw("l");
  TLatex *text = new TLatex();
  text->SetTextSize(0.05);
  text->SetTextColor(4);
  text->DrawLatex(1,42,"Electron");
  text->SetTextColor(2);
  text->DrawLatex(40,30,"Muon");
  c1->SaveAs("hist_PhD/opening_angle.pdf");
  /*TLegend *leg = new TLegend(0.1,0.6,0.3,0.9);
  leg->SetFillColor(0);
  leg->SetBorderSize(0);
  leg->AddEntry(g_e,"Electron","l");
  leg->AddEntry(g_mu,"Muon","l");
  leg->Draw();*/


}

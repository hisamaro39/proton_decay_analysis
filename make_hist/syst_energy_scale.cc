#include <vector>
void syst_energy_scale(){
  string mode = "p_eee";
  int cut = 5;
  int nring=1;
  int mulike=0;
  int michel=0;

  gROOT->SetStyle("Plain");
  gStyle->SetOptStat(0);

  cout << "mode=" << mode << endl;
  TH1 *first_hist;
  TFile *input = TFile::Open(Form("../output/fcmc.sk4.mode_%s.root",mode.c_str()));//mc
  TH1* h_total_distance = (TH1*) input->Get(Form("total_distance_ntag_cut%d_nring%d_mulike%d_michel%d",cut,nring,mulike,michel));
  TCanvas *c = new TCanvas("c","",800,600);
  h_total_distance->SetLineWidth(2);
  h_total_distance->Draw("hist");
  TF1 *func = new TF1("func","exp([0]*x*x+[1]*x+[2])",0,350);
  //TF1 *func = new TF1("func","[0]*x*x*x+[1]*x*x+[2]*x",0,500);
  func->SetLineColor(2);
  h_total_distance->Fit(func,"","",0,350);
  float par0 = func->GetParameter(0);
  float par1 = func->GetParameter(1);
  float par2 = func->GetParameter(2);
  cout << "function is " << "exp(" << par0 << "*x^2 + " << par1 << "*x + " << par2 << ")" << endl;
  TF1 *fix_func = new TF1("fix_func","exp([0]*x*x+[1]*x+[2])",0,350);
  fix_func->FixParameter(0,par0);
  fix_func->FixParameter(1,par1);
  fix_func->FixParameter(2,par2);
  fix_func->SetLineColor(4);
  fix_func->Draw("same");
  c->SaveAs(Form("hist/total_distance_ntag_cut%d_nring%d_mulike%d_michel%d_fit.pdf",cut,nring,mulike,michel));
  double total_events = fix_func->Integral(0,125)/25;
  double total_events_up = fix_func->Integral(0,125*(1+0.021))/25;
  double total_events_down = fix_func->Integral(0,125*(1-0.021))/25;
  double diff_up = (total_events_up - total_events)/total_events;
  double diff_down = (total_events_down - total_events)/total_events;
  float err2 = 0.;
  for (int b=1;b<=5;b++) err2 += pow(h_total_distance->GetBinError(b),2);
  float err = sqrt(err2);
  cout << "integral=" << h_total_distance->Integral(1,5) << " +- " << err << endl;
  cout << "total events def/up/down=" << total_events << "/" << total_events_up << "/" << total_events_down << endl;
  cout << "diff up/down=" << diff_up << "/" << diff_down << endl; 

  /*TF1 *test_func = new TF1("test_func","[0]*x",0,600);
  test_func->FixParameter(0,1);
  test_func->SetLineColor(3);
  test_func->Draw();
  cout << test_func->Integral(0,10) << endl;*/
}

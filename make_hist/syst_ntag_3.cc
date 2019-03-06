#include <vector>
void syst_ntag_3(){
  string mode = "p_mumumu";
  int iteration = 50;
  int cut = 7;
  int nring=0;
  int mulike=0;
  int michel=3;

  gROOT->SetStyle("Plain");
  gStyle->SetOptStat(0);

  cout << "mode=" << mode << endl;
  TH1 *first_hist;
  TFile *input_mc = TFile::Open(Form("../output/fcmc.sk4.mode_%s.root",mode.c_str()));//mc
  TH1* h_n_neutron_true = (TH1*) input_mc->Get(Form("n_true_neutron_cut5_nring%d_mulike%d_michel%d",nring,mulike,michel));
  TH1* h_n_neutron_tag = (TH1*) input_mc->Get(Form("ntag_multiplicity_cut%d_nring%d_mulike%d_michel%d",cut,nring,mulike,michel));
  //TCanvas *c1 = new TCanvas("c1","",800,600);
  h_n_neutron_tag->SetLineWidth(2);
  h_n_neutron_true->SetLineWidth(2);
  h_n_neutron_true->SetLineColor(2);
  //h_n_neutron_tag->Draw();
  //h_n_neutron_true->Draw("same");
  //c1->SaveAs("hist/compare_true_neutron_tagged_neutron_cut5_nring1_mulike0_michel0_p_eee_sk4.pdf");
  float n_no_tagged = h_n_neutron_tag->GetBinContent(1);
  float n_no_tagged_calc=0.;
  float efficiency = 0.2;
  TRandom *generator = new TRandom();
  generator->SetSeed(0);
  for(int b=0;b<10;b++){
    int n_neutron_true = h_n_neutron_true->GetBinContent(b+1);
    float n_no_tagged_this = pow(1-efficiency,b) * n_neutron_true;
    //cout << "n_no_tagged_this=" << n_no_tagged_this << endl;
    n_no_tagged_calc += n_no_tagged_this;
    float n_events_this=0.;
    //for(int bb=0;bb<b;bb++){
      //float cont_this = pow(1-efficiency,bb-b)*pow(efficiency,)
    //}
  }
  //cout << "no ntag events detect/calc=" << n_no_tagged << "/" << n_no_tagged_calc << endl;

  TFile *input = TFile::Open(Form("../output/fcmc.sk4.mode_%s_syst_ntag_it%d.root",mode.c_str(),iteration));//mc
  TH1* h_n_neutron_true_final = (TH1*) input_mc->Get(Form("n_true_neutron_cut%d_nring%d_mulike%d_michel%d",cut,nring,mulike,michel));
  float total_events = h_n_neutron_true_final->Integral();
  float total_err2 = 0.;
  for(int b=1;b<=h_n_neutron_true_final->GetNbinsX();b++)
    total_err2 += pow(h_n_neutron_true_final->GetBinError(b),2);
  float total_err = sqrt(total_err2)/total_events;
  cout << "total events = " << total_events << " += " << sqrt(total_err2) << endl;
  cout << "MC stat error is " << total_err << endl;

  TGraphErrors *graph = new TGraphErrors();
  for(int e=0;e<11;e++){
    int eff = 10*e;
    cout << "assumed efficiency is " << eff << "%" << endl;
    TH1* h_n_neutron_exp = (TH1*) input->Get(Form("n_tagged_neutron_exp_eff%d_cut%d_nring%d_mulike%d_michel%d",eff,cut,nring,mulike,michel));
    float total_events = h_n_neutron_exp->Integral();
    float sum_2err = 0;
    for(int b=0;b<h_n_neutron_exp->GetNbinsX();b++) sum_2err += pow(h_n_neutron_exp->GetBinError(b+1),2);
    float err_total_events = sqrt(sum_2err);
    float no_ntag_events = h_n_neutron_exp->GetBinContent(1);
    float err_no_ntag_events = h_n_neutron_exp->GetBinError(1);
    float fraction = no_ntag_events / total_events;
    float err = sqrt(pow(no_ntag_events*err_total_events,2)+pow(total_events*err_no_ntag_events,2))/pow(total_events,2);
    cout << "total_events=" << total_events << endl;
    cout << "no_ntag_events=" << no_ntag_events << endl;
    cout << "fraction is " << fraction << " +- " << fraction*total_err << endl;
    graph->SetPoint(e,0.1*e,fraction);
    graph->SetPointError(e,0,fraction*total_err);
  }
  TCanvas *c2 = new TCanvas("c2","",800,600);
  graph->SetMarkerStyle(8);
  TH1 *frame=gPad->DrawFrame(-0.1, 0, 1.1, 1.1);
  graph->Draw("p");
  TF1 *func = new TF1("func","[0]+exp([1]*x+[2])",0,1);
  func->SetLineColor(2);
  graph->Fit(func);
  c2->SaveAs(Form("hist/ntag_eff_bkg_frac_cut7_%s_sk4.pdf",mode.c_str()));
  float par0 = func->GetParameter(0);
  float par1 = func->GetParameter(1);
  float par2 = func->GetParameter(2);
  cout << "fit result = " << par0 << " + exp(" << par1 << "x + " << par2 << ")" << endl;
  float frac_nom = par0+exp(par1*0.205+par2);
  float frac_down = par0+exp(par1*0.1845+par2);
  float frac_up = par0+exp(par1*0.2255+par2);
  cout << "bkg frac for eff 20.5% is " << frac_nom << endl;
  cout << "bkg frac for eff 18.45% is " << frac_down << endl;
  cout << "bkg frac for eff 22.55% is " << frac_up << endl;
  float diff_up = (frac_nom - frac_up)/frac_nom; 
  float diff_down = (frac_nom - frac_down)/frac_nom; 
  cout << "difference up/down=" << diff_up << "/" << diff_down << endl;

  TFile *input_it10 = TFile::Open(Form("../output/fcmc.sk4.mode_%s_syst_ntag_it10.root",mode.c_str()));//mc
  TFile *input_it30 = TFile::Open(Form("../output/fcmc.sk4.mode_%s_syst_ntag_it30.root",mode.c_str()));//mc
  TFile *input_it50 = TFile::Open(Form("../output/fcmc.sk4.mode_%s_syst_ntag_it50.root",mode.c_str()));//mc
  TCanvas *c3 = new TCanvas("c3","",800,600);
  TH1* hist_eff20_it10 = (TH1*) input_it10->Get(Form("n_tagged_neutron_exp_eff20_cut%d_nring%d_mulike%d_michel%d",cut,nring,mulike,michel));
  TH1* hist_eff20_it30 = (TH1*) input_it30->Get(Form("n_tagged_neutron_exp_eff20_cut%d_nring%d_mulike%d_michel%d",cut,nring,mulike,michel));
  TH1* hist_eff20_it50 = (TH1*) input_it50->Get(Form("n_tagged_neutron_exp_eff20_cut%d_nring%d_mulike%d_michel%d",cut,nring,mulike,michel));
  hist_eff20_it10->SetLineWidth(2);
  hist_eff20_it30->SetLineWidth(2);
  hist_eff20_it50->SetLineWidth(2);
  hist_eff20_it10->SetLineColor(2);
  hist_eff20_it30->SetLineColor(3);
  hist_eff20_it50->SetLineColor(2);
  h_n_neutron_tag->Draw("hist E0");
  hist_eff20_it10->Scale(1./10);
  hist_eff20_it30->Scale(1./30);
  hist_eff20_it50->Scale(1./50);
  //hist_eff20_it10->Draw("hist same");
  //hist_eff20_it30->Draw("hist same");
  hist_eff20_it50->Draw("hist same");
  c3->SaveAs(Form("hist/compare_ntag_reco_20_cut%d_%s_sk4.pdf",cut,mode.c_str()));

}

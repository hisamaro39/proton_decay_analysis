#include <vector>
void syst_ntag_3(){
  string mode = "p_eee";

  gROOT->SetStyle("Plain");
  gStyle->SetOptStat(0);

  cout << "mode=" << mode << endl;
  TH1 *first_hist;
  TFile *input_mc = TFile::Open(Form("../output/fcmc.sk4.mode_%s.root",mode.c_str()));//mc
  TH1* h_n_neutron_true = (TH1*) input_mc->Get("n_true_neutron_cut5_nring1_mulike0_michel0");
  TH1* h_n_neutron_tag = (TH1*) input_mc->Get("ntag_multiplicity_cut5_nring1_mulike0_michel0");
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

  TFile *input = TFile::Open(Form("../output/fcmc.sk4.mode_%s_syst_ntag.root",mode.c_str()));//mc
  TH1* h_n_neutron_true_final = (TH1*) input_mc->Get("n_true_neutron_cut7_nring1_mulike0_michel0");
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
    TH1* h_n_neutron_exp = (TH1*) input->Get(Form("n_tagged_neutron_exp_eff%d_cut7_nring1_mulike0_michel0",eff));
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
  c2->SaveAs("hist/ntag_eff_bkg_frac_cut7_p_eee_sk4.pdf");
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

  //TCanvas *c3 = new TCanvas("c3","",800,600);
  TH1* hist_10 = (TH1*) input->Get("n_tagged_neutron_exp_eff10_cut5_nring1_mulike0_michel0");
  TH1* hist_20 = (TH1*) input->Get("n_tagged_neutron_exp_eff20_cut5_nring1_mulike0_michel0");
  TH1* hist_30 = (TH1*) input->Get("n_tagged_neutron_exp_eff30_cut5_nring1_mulike0_michel0");
  hist_10->SetLineWidth(2);
  hist_20->SetLineWidth(2);
  hist_30->SetLineWidth(2);
  hist_10->SetLineColor(2);
  hist_20->SetLineColor(3);
  hist_30->SetLineColor(4);
  //h_n_neutron_tag->GetYaxis()->SetRangeUser(0,400);
  //h_n_neutron_tag->Draw();
  //hist_10->Draw("same");
  //hist_20->Draw("same");
  //hist_30->Draw("same");
  //c3->SaveAs("hist/compare_exp_ntag_real_10_20_30_cut5_p_eee_sk4.pdf");*/

  /*int no_cap_n = hist_low->GetBinContent(1) + hist_high->GetBinContent(1);
  int cap_n = hist_low->GetBinContent(2) + hist_high->GetBinContent(2);
  cout << "# of events no_captured/captured=" << no_cap_n << "/" << cap_n << endl;

  TGraphErrors *graph = new TGraphErrors();
  TRandom *generator = new TRandom();
  int n_iteration = 10;
  for(int p=0;p<11;p++){
    float ntag_eff = 0.1*p;
    cout << "assumed ntag efficiency = " << ntag_eff << endl;
    int tagged_cap_n=0;
    for(int i=0;i<n_iteration;i++){
      cout << "iteration:" << i << endl;
      generator->SetSeed(i);
      for(int n=0;n<cap_n;n++){
        float rnd = generator->Rndm();
        cout << "rnd=" << rnd << endl;
        if(rnd<ntag_eff) {
          tagged_cap_n++;
          cout << "tag!!" << endl;
        }
      }
      cout << "tagged_cap_n=" << tagged_cap_n << endl;
    }
    int total_cap_n = cap_n * n_iteration;
    cout << "cap_n/tagged=" << total_cap_n << "/" << tagged_cap_n << endl;
    float frac = 1.*tagged_cap_n / total_cap_n;
    float frac_err = (tagged_cap_n)? frac*sqrt(1./total_cap_n + 1./tagged_cap_n) : frac/sqrt(total_cap_n);
    cout << "efficiency is " << frac << " +- " << frac_err << endl;
    int total_initial = n_iteration * ( no_cap_n + cap_n );
    int total_remained = total_initial - tagged_cap_n;
    float total_frac = 1.*total_remained/total_initial;
    float total_frac_err = total_frac*sqrt(1./total_remained + 1./total_initial);
    cout << "total fraction is " << total_frac << " +- " << total_frac_err << endl;
    graph->SetPoint(p,ntag_eff,total_frac);
    graph->SetPointError(p,0,total_frac_err);
  }
  graph->SetMarkerStyle(8);
  graph->Draw("ap");
  */

}


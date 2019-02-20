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
    cout << "n_no_tagged_this=" << n_no_tagged_this << endl;
    n_no_tagged_calc += n_no_tagged_this;
    float n_events_this=0.;
    //for(int bb=0;bb<b;bb++){
      //float cont_this = pow(1-efficiency,bb-b)*pow(efficiency,)
    //}
  }
  cout << "no ntag events detect/calc=" << n_no_tagged << "/" << n_no_tagged_calc << endl;

  TFile *input = TFile::Open(Form("../output/fcmc.sk4.mode_%s_syst_ntag.root",mode.c_str()));//mc
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
    cout << "fraction is " << fraction << " +- " << err << endl;
    graph->SetPoint(e+1,0.1*e,fraction);
    graph->SetPointError(e+1,0,err);
  }
  graph->SetMarkerStyle(8);
  graph->Draw("ap");

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


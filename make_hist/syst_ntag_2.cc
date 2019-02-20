#include <vector>
void syst_ntag_2(){
  string mode = "p_eee";

  gROOT->SetStyle("Plain");
  gStyle->SetOptStat(0);

  cout << "mode=" << mode << endl;
  TH1 *first_hist;
  TFile *input_mc = TFile::Open(Form("../output/fcmc.sk4.mode_%s.root",mode.c_str()));//mc
  TH1* hist_low = (TH1*) input_mc->Get("has_captured_neutron_cut6_nring1_mulike0_michel0");//low SR w/o ntag
  TH1* hist_high = (TH1*) input_mc->Get("has_captured_neutron_cut7_nring1_mulike0_michel0");//high SR w/o ntag
  int no_cap_n = hist_low->GetBinContent(1) + hist_high->GetBinContent(1);
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

}


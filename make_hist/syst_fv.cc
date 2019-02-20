#include <vector>
void syst_fv(){
  string type[] = {"p_epi","p_mupi","p_eee","p_mumumu","p_emumu","p_muee","fcmc","fcdt","subgev_multiring"};
  int mode_id = 8;
  string sk_period = "sk4";
  bool add_signal = false;
  float signal_scale = 0.03;

  int range_momentum[] = {100,250,400,630,1000,2500,5000,10000,100000};
  int n_range_momentum = sizeof(range_momentum)/sizeof(int);

  gROOT->SetStyle("Plain");
  gStyle->SetOptStat(0);
  vector<string> hist;
  vector<vector<string>> hist_set;
  vector<vector<int>> input_type_set,mode_type_set;
  vector<int> scale,dology,input_type,mode_type,add_ratio,dorebin;
  hist.clear();hist_set.clear();scale.clear();dology.clear();dorebin.clear();input_type.clear();input_type_set.clear();add_ratio.clear();
  mode_type_set.clear();mode_type.clear();
  vector<string> hist_name;
  hist_name.clear();

  TH1 *first_hist;
  TFile *input_mc = TFile::Open("../output/fcmc.sk4.mode_subgev_multiring.root");//mc
  TFile *input_data = TFile::Open("../output/fcdt.sk4.mode_subgev_multiring.root");//mc
  int dwall_thr[3] = {50,100,150};
  for(int d=0;d<3;d++){
    cout << "dwall_thr=" << dwall_thr[d] << endl;
    TH1* hist_mc = (TH1*) input_mc->Get(Form("distance_to_wall_thr%d",dwall_thr[d]));
    TH1* hist_data = (TH1*) input_data->Get(Form("distance_to_wall_thr%d",dwall_thr[d]));
    hist_data->SetLineWidth(2);
    hist_mc->SetLineWidth(2);
    //cout << hist_mc->GetEntries() << "/" << hist_data->GetEntries() << endl;
    TCanvas *c1 = new TCanvas(Form("c1_%d",d),"",800,600);
    hist_data->Draw();
    hist_mc->SetLineColor(2);
    hist_mc->Draw("same hist");
    c1->SaveAs(Form("hist/compare_data_mc_distance_to_wall_thr%d.pdf",dwall_thr[d]));
    float evt_dwall_thr_data = hist_data->Integral(d+2,36);
    float evt_dwall_thr_mc = hist_mc->Integral(d+2,36);
    //cout << "events in dwall_thr(integral) data/mc=" << evt_dwall_thr_data << "/" << evt_dwall_thr_mc << endl;
    float norm_factor = evt_dwall_thr_data/evt_dwall_thr_mc;
    //cout << "norm_factor=" << norm_factor << endl;
    TH1* norm_hist_mc = (TH1*) hist_mc->Clone("norm_hist_mc");
    norm_hist_mc->Scale(norm_factor);
    norm_hist_mc->SetLineColor(4);
    norm_hist_mc->SetLineWidth(2);
    TCanvas *c2 = new TCanvas(Form("c2_%d",d),"",800,600);
    hist_data->Draw();
    norm_hist_mc->Draw("same hist");
    c2->SaveAs(Form("hist/compare_data_mc_distance_to_wall_scaled_thr%d.pdf",dwall_thr[d]));
    float evt_dwall_thr_scale_mc = norm_hist_mc->Integral(d+2,36);
    float evt_fv_data = hist_data->Integral(5,36);
    float evt_fv_mc = norm_hist_mc->Integral(5,36);
    //cout << "events in dwall_thr data/scale_mc=" << evt_dwall_thr_data << "/" << evt_dwall_thr_scale_mc << endl;
    //cout << "events in FV data/scale_mc=" << evt_fv_data << "/" << evt_fv_mc << endl;
    float sum_mc=0.,sum_data=0.,sum_mc_norm=0.,sumerr_mc=0.,sumerr_data=0.,sumerr_mc_norm=0.;;
    for(int b=5;b<37;b++){
      float evt_mc = hist_mc->GetBinContent(b);
      float err_mc = hist_mc->GetBinError(b);
      float evt_mc_norm = norm_hist_mc->GetBinContent(b);
      float err_mc_norm = norm_hist_mc->GetBinError(b);
      float evt_data = hist_data->GetBinContent(b);
      float err_data = hist_data->GetBinError(b);
      //cout << "bin" << b << " evt_mc=" << evt_mc << " +- " << err_mc << endl;
      //cout << "bin" << b << " evt_mc_norm=" << evt_mc_norm << " +- " << err_mc_norm << endl;
      //cout << "bin" << b << " evt_data=" << evt_data << " +- " << err_data << endl;
      sum_mc += evt_mc;
      sum_mc_norm += evt_mc_norm;
      sumerr_mc += err_mc*err_mc;
      sumerr_mc_norm += err_mc_norm*err_mc_norm;
      sum_data += evt_data;
      sumerr_data += err_data*err_data;
    }
    cout << "events in FV data =" << sum_data << endl;
    //cout << "events in FV mc =" << sum_mc << " +- " << sqrt(sumerr_mc) << endl;
    cout << "events in FV mc_norm =" << sum_mc_norm << " +- " << sqrt(sumerr_mc_norm) << endl;
    float diff = sum_data - sum_mc_norm;
    float diff_err = sqrt(sumerr_mc_norm);
    float diff_ratio = diff / sum_data;
    float diff_ratio_err = sqrt(sumerr_mc_norm) / sum_data;
    cout << "data - mc = " << diff << " +- " << diff_err << endl;
    cout << "diff / data = " << diff_ratio << " +- " << diff_ratio_err << endl; 
  }
}

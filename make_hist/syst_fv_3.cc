#include <vector>
void syst_fv_3(){
  string type[] = {"p_epi","p_mupi","p_eee","p_mumumu","p_emumu","p_muee","fcmc","fcdt","subgev_multiring"};
  int mode_id = 8;
  string sk_period = "sk4";
  bool add_signal = false;
  float signal_scale = 0.03;
  //string hist_name = "distance_to_wall_thr50";
  string hist_name = "vertex_z";

  gROOT->SetStyle("Plain");
  gStyle->SetOptStat(0);
  vector<string> hist;
  vector<vector<string>> hist_set;
  vector<vector<int>> input_type_set,mode_type_set;
  vector<int> scale,dology,input_type,mode_type,add_ratio,dorebin;
  hist.clear();hist_set.clear();scale.clear();dology.clear();dorebin.clear();input_type.clear();input_type_set.clear();add_ratio.clear();
  mode_type_set.clear();mode_type.clear();

  TH1 *first_hist;
  TFile *input_mc_sk2 = TFile::Open("../output/fcmc.sk2.mode_subgev_multiring.root");//mc
  TFile *input_mc_sk4 = TFile::Open("../output/fcmc.sk4.mode_subgev_multiring.root");//mc
  TFile *input_data_sk2 = TFile::Open("../output/fcdt.sk2.mode_subgev_multiring.root");//data
  TFile *input_data_sk4 = TFile::Open("../output/fcdt.sk4.mode_subgev_multiring.root");//data
  
  //MC
  TH1* norm_hist_mc_sk2 = (TH1*) input_mc_sk2->Get("distance_to_wall_thr50");
  TH1* norm_hist_mc_sk4 = (TH1*) input_mc_sk4->Get("distance_to_wall_thr50");
  float evt_mc_sk2 = norm_hist_mc_sk2->Integral(13,20);//600 < Dwall < 1000 cm
  float evt_mc_sk4 = norm_hist_mc_sk4->Integral(13,20);//600 < Dwall < 1000 cm
  cout << "events in 600<Dwall<1200cm mc sk2/sk4=" << evt_mc_sk2 << "/" << evt_mc_sk4 << endl;
  float norm_factor_mc = evt_mc_sk4/evt_mc_sk2;
  TH1* hist_mc_sk2 = (TH1*) input_mc_sk2->Get(hist_name.c_str());
  TH1* hist_mc_sk4 = (TH1*) input_mc_sk4->Get(hist_name.c_str());
  hist_mc_sk2->Scale(norm_factor_mc);
  TCanvas *c1 = new TCanvas("c1","",800,600);
  TPad* p1 = new TPad("main","main",0.0,0.2,1.0,1.0,10,0,0);
  p1->SetNumber(1);
  p1->SetBottomMargin(0);
  p1->Draw();
  TPad* p2 = new TPad("ratio","ratio",0.0,0.0,1.0,0.2,10,0,0);
  p2->SetTopMargin(0);
  p2->SetBottomMargin(0.5);
  p2->SetNumber(2);
  p2->Draw();
  c1->cd(1);
  hist_mc_sk2->SetLineWidth(2);
  hist_mc_sk4->SetLineWidth(2);
  hist_mc_sk4->SetMinimum(0);
  hist_mc_sk4->Draw("hist E0");
  hist_mc_sk2->SetLineColor(2);
  hist_mc_sk2->Draw("same hist E0");
  TH1 *ratio_hist_mc = (TH1*) hist_mc_sk2->Clone("clone_hist_mc_sk2");
  ratio_hist_mc->Divide(hist_mc_sk4);
  float xmin = ratio_hist_mc->GetBinLowEdge(1);
  float xmax = ratio_hist_mc->GetBinLowEdge(ratio_hist_mc->GetNbinsX())+ratio_hist_mc->GetBinWidth(ratio_hist_mc->GetNbinsX());
  c1->cd(2);
  TH1* frame;
  frame=gPad->DrawFrame(xmin, 0.5, xmax, 1.5);
  frame->GetYaxis()->SetLabelSize(0.1);
  frame->GetXaxis()->SetLabelSize(0.2);
  ratio_hist_mc->Draw("same");
  TLine *line = new TLine(xmin,1,xmax,1);
  line->SetLineStyle(2);
  line->SetLineWidth(2);
  line->Draw();
  c1->SaveAs(Form("hist/compare_%s_mc_sk2_sk4.pdf",hist_name.c_str()));

  //Data
  TH1* norm_hist_data_sk2 = (TH1*) input_data_sk2->Get("distance_to_wall_thr50");
  TH1* norm_hist_data_sk4 = (TH1*) input_data_sk4->Get("distance_to_wall_thr50");
  float evt_data_sk2 = norm_hist_data_sk2->Integral(13,20);//600 < Dwall < 1000 cm
  float evt_data_sk4 = norm_hist_data_sk4->Integral(13,20);//600 < Dwall < 1000 cm
  cout << "events in 600<Dwall<1200cm data sk2/sk4=" << evt_data_sk2 << "/" << evt_data_sk4 << endl;
  float norm_factor_data = evt_data_sk4/evt_data_sk2;
  TH1* hist_data_sk2 = (TH1*) input_data_sk2->Get(hist_name.c_str());
  TH1* hist_data_sk4 = (TH1*) input_data_sk4->Get(hist_name.c_str());
  hist_data_sk2->Scale(norm_factor_data);
  TCanvas *c2 = new TCanvas("c2","",800,600);
  TPad* p3 = new TPad("main","main",0.0,0.2,1.0,1.0,10,0,0);
  p3->SetNumber(1);
  p3->SetBottomMargin(0);
  p3->Draw();
  TPad* p4 = new TPad("ratio","ratio",0.0,0.0,1.0,0.2,10,0,0);
  p4->SetTopMargin(0);
  p4->SetBottomMargin(0.5);
  p4->SetNumber(2);
  p4->Draw();
  c2->cd(1);
  float evt_max = (hist_data_sk2->GetMaximum() > hist_data_sk4->GetMaximum())?
    hist_data_sk2->GetMaximum() : hist_data_sk4->GetMaximum();
  hist_data_sk4->GetYaxis()->SetRangeUser(0,evt_max*1.2);
  hist_data_sk2->SetLineWidth(2);
  hist_data_sk4->SetLineWidth(2);
  hist_data_sk4->Draw("hist E0");
  hist_data_sk2->SetLineColor(2);
  hist_data_sk2->Draw("same hist E0");
  TH1 *ratio_hist_data = (TH1*) hist_data_sk2->Clone("clone_hist_data_sk2");
  ratio_hist_data->Divide(hist_data_sk4);
  c2->cd(2);
  TH1* frame2;
  float xmin2 = ratio_hist_data->GetBinLowEdge(1);
  float xmax2 = ratio_hist_data->GetBinLowEdge(ratio_hist_data->GetNbinsX())+ratio_hist_data->GetBinWidth(ratio_hist_data->GetNbinsX());
  frame2=gPad->DrawFrame(xmin2, 0.5, xmax2, 1.5);
  frame2->GetYaxis()->SetLabelSize(0.1);
  frame2->GetXaxis()->SetLabelSize(0.2);
  ratio_hist_data->Draw("same");
  TLine *line2 = new TLine(xmin2,1,xmax2,1);
  line2->SetLineStyle(2);
  line2->SetLineWidth(2);
  line2->Draw();
  c2->SaveAs(Form("hist/compare_%s_data_sk2_sk4.pdf",hist_name.c_str()));
}

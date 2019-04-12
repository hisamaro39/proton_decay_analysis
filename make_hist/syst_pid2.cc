void syst_pid2(){

  gROOT->SetStyle("Plain");
  gStyle->SetOptStat(0);

  TFile *input_mc = TFile::Open("../output/fcmc.sk4.mode_subgev_multiring.root");//mc
  TFile *input_data = TFile::Open("../output/fcdt.sk4.mode_subgev_multiring.root");//mc
  TH1* hist_mc = (TH1*) input_mc->Get("prob_angle_3rd_nring3");
  TH1* hist_data = (TH1*) input_data->Get("prob_angle_3rd_nring3");

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
  hist_data->SetLineWidth(2);
  hist_mc->SetLineWidth(2);
  hist_mc->SetLineColor(2);
  hist_data->Draw();
  hist_mc->Draw("same hist");
  TH1 *ratio_hist_data = (TH1*) hist_data->Clone("clone_hist_data");
  ratio_hist_data->Divide(hist_mc);
  float xmin = ratio_hist_data->GetBinLowEdge(1);
  float xmax = ratio_hist_data->GetBinLowEdge(ratio_hist_data->GetNbinsX())+ratio_hist_data->GetBinWidth(ratio_hist_data->GetNbinsX());
  c1->cd(2);
  TH1* frame;
  frame=gPad->DrawFrame(xmin, 0.5, xmax, 1.5);
  frame->GetYaxis()->SetLabelSize(0.1);
  frame->GetXaxis()->SetLabelSize(0.2);
  ratio_hist_data->Draw("same");
  TLine *line = new TLine(xmin,1,xmax,1);
  line->SetLineStyle(2);
  line->SetLineWidth(2);
  line->Draw();
  c1->SaveAs("hist/compare_data_mc_prob_angle_subgev_multiring_3rd.pdf");

  //calculate chi2
  int nbin = hist_data->GetNbinsX();
  cout << "# of bins " << nbin << endl;
  TGraph *graph = new TGraph();
  TGraph *graph2 = new TGraph();
  float total_min_chi2 = 999999,total_min_scale=0;
  int total_min_shift=9999;
  for(int m=0;m<10;m++){//bin shift
    int shift=-5+m; 
    cout << "bin shift is " << shift << endl;
    float min_chi2 = 999999,min_scale=0,min_shift=9999;
    for(int s=0;s<1000;s++){//scale
      float scale = 0.5 + 0.001*s;
      //cout << "scale is " << scale << endl;
      float chi2=0.;
      for(int b=0;b<nbin;b++){
        float event_data = hist_data->GetBinContent(b+1);
        float event_mc = hist_mc->GetBinContent(b+1+shift)*scale;
        float sigma = sqrt(event_data);
        if(event_data==0) continue;
        float this_diff = pow((event_data-event_mc)/sigma,2);
        chi2 += this_diff;
        //cout << "event_data/event_mc/sigma=" << event_data << "/" << event_mc << "/" << sigma << endl;
      }
      chi2 = chi2/nbin;
      if(shift==0)graph->SetPoint(s,scale,chi2);//make example plot
      if(chi2<min_chi2) {
        min_chi2 = chi2;
        min_scale = scale;
        min_shift = shift;
      }
      //graph->SetPoint(s,scale,chi2);
      //cout << "chi2=" << chi2 << endl;
    }
    cout << "minimum chi2/shift/scale " << min_chi2 << "/" << min_shift << "/" << min_scale << endl;
    if(min_chi2<total_min_chi2){
      total_min_chi2 = min_chi2;
      total_min_shift = min_shift;
      total_min_scale = min_scale;
    }
    graph2->SetPoint(m,shift,min_chi2);
  }
  cout << "total minimum chi2/shift/scale=" << total_min_chi2 << "/" << total_min_shift << "/" << total_min_scale << endl;

  TCanvas *c4 = new TCanvas("c4","",800,600);
  graph->SetLineWidth(2);
  graph->Draw("al");
  c4->SaveAs("hist/chi2_shift0_prob_angle_subgev_multiring_3rd.pdf");

  TCanvas *c2 = new TCanvas("c2","",800,600);
  TF1 *func = new TF1("func","[0]*x*x+[1]*x+[2]",-3,3);
  func->SetLineColor(2);
  graph2->Fit(func,"","",-2,2);
  graph2->SetMarkerStyle(8);
  graph2->Draw("ap");
  c2->SaveAs("hist/minchi2_prob_angle_subgev_multiring_3rd.pdf");
  float par0 = func->GetParameter(0);
  float par1 = func->GetParameter(1);
  float par2 = func->GetParameter(2);
  cout << "fitted function " << par0 << "*x^2+" << par1 << "*x+" << par2 << endl;
  TF1 *fix_func = new TF1("fix_func","[0]*x*x+[1]*x+[2]",-2,2);
  fix_func->FixParameter(0,par0);
  fix_func->FixParameter(1,par1);
  fix_func->FixParameter(2,par2);
  fix_func->SetLineColor(2);
  fix_func->SetLineWidth(2);
  //fix_func->Draw("same");
  float min_value=9999,min_x=9999;
  for(int i=0;i<4000;i++){
    float x = -2+0.001*i;
    float value = fix_func->Eval(x);
    if(value<min_value) {
      min_value = value;
      min_x = x;
    }
  }
  cout << "minimum chi2/bin=" << min_value << "/" << min_x << endl;
  float one_sigma_bin = 9999;
  for(int j=0;j<3000;j++){
    float x2 = -3+0.001*j;
    float value2 = fix_func->Eval(x2);
    //cout << "x2/value2-min=" << x2 << "/" << value2-min_value << endl;
    if(value2-min_value<1) {
      one_sigma_bin=x2;
      break;
    }
  }
  float one_sigma_bin2 = 9999;
  for(int k=0;k<3000;k++){
    float x3 = 0.001*k;
    float value3 = fix_func->Eval(x3);
    //cout << "x2/value2-min=" << x2 << "/" << value2-min_value << endl;
    if(value3-min_value>1) {
      one_sigma_bin2=x3;
      break;
    }
  }
  cout << "one sigma bin is " << one_sigma_bin << "/" << one_sigma_bin2 << endl;
  cout << "diff from minimum is " << one_sigma_bin - min_x << "/" << one_sigma_bin2 - min_x << endl;
  float bin_width = hist_data->GetBinWidth(1);
  cout << "prob_angle shift value for 1 sigma is " << bin_width*(one_sigma_bin-min_x) 
    << "/" << bin_width*(one_sigma_bin2-min_x) << endl;

  TCanvas *c3 = new TCanvas("c3","",800,600);
  TPad* p3 = new TPad("main","main",0.0,0.2,1.0,1.0,10,0,0);
  p3->SetNumber(1);
  p3->SetBottomMargin(0);
  p3->Draw();
  TPad* p4 = new TPad("ratio","ratio",0.0,0.0,1.0,0.2,10,0,0);
  p4->SetTopMargin(0);
  p4->SetBottomMargin(0.5);
  p4->SetNumber(2);
  p4->Draw();
  c3->cd(1);
  TH1 *hist_mc_fit = (TH1*) hist_mc->Clone("clone_hist_mc");
  hist_mc_fit->Scale(0);
  hist_mc_fit->Draw();
  for(int b=0;b<nbin;b++){
    float shifted_value = hist_mc->GetBinContent(b+1+total_min_shift);
    hist_mc_fit->SetBinContent(b+1,shifted_value);
  }
  hist_mc_fit->Scale(total_min_scale);
  hist_data->SetLineWidth(2);
  hist_mc_fit->SetLineWidth(2);
  hist_mc_fit->SetLineColor(4);
  hist_data->Draw();
  hist_mc_fit->Draw("same hist");
  TH1 *ratio_hist_data2 = (TH1*) hist_data->Clone("clone_hist_datai2");
  ratio_hist_data2->Divide(hist_mc_fit);
  c3->cd(2);
  TH1* frame2;
  frame2=gPad->DrawFrame(xmin, 0.5, xmax, 1.5);
  frame2->GetYaxis()->SetLabelSize(0.1);
  frame2->GetXaxis()->SetLabelSize(0.2);
  ratio_hist_data2->Draw("same");
  TLine *line2 = new TLine(xmin,1,xmax,1);
  line2->SetLineStyle(2);
  line2->SetLineWidth(2);
  line2->Draw();
  c3->SaveAs("hist/compare_data_mc_prob_angle_bestfit_subgev_multiring_3rd.pdf");

  TCanvas *c5 = new TCanvas("c5","",800,600);
  hist_data->Draw();
  hist_mc->Draw("hist same");
  hist_mc_fit->Draw("hist same");
  c5->SaveAs("hist/compare_data_mc_prob_angle_before_after_fit_subgev_multiring_3rd.pdf");

}

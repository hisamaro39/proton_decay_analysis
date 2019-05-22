#include <vector>
void syst_energy_scale(){
  string mode = "p_mumumu";
  int cut = 5;
  int nring=0;
  int mulike=0;
  int michel=3;
  int fit_max=600;

  gROOT->SetStyle("Plain");
  gStyle->SetOptStat(0);

  cout << "mode=" << mode << endl;
  TH1 *first_hist;
  TFile *input = TFile::Open(Form("../output/fcmc_real.sk1_4.mode_%s_wo_livetime.root",mode.c_str()));//mc
  TH1* h_total_distance = (TH1*) input->Get(Form("total_distance_ntag_cut%d_nring%d_mulike%d_michel%d",cut,nring,mulike,michel));
  TCanvas *c1 = new TCanvas("c1","",800,600);
  h_total_distance->SetLineWidth(2);
  h_total_distance->Draw("hist");
  TF1 *func = new TF1("func","exp([0]*x*x+[1]*x+[2])",0,fit_max);
  func->SetLineColor(2);
  h_total_distance->Fit(func,"","",0,fit_max);
  float par0 = func->GetParameter(0);
  float par1 = func->GetParameter(1);
  float par2 = func->GetParameter(2);
  float par0_err = func->GetParError(0);
  float par1_err = func->GetParError(1);
  float par2_err = func->GetParError(2);
  float pm_par0_err[2] = {par0_err,-1*par0_err};
  float pm_par1_err[2] = {par1_err,-1*par1_err};
  float pm_par2_err[2] = {par2_err,-1*par2_err};
  cout << "function is " << "exp(" << par0 << "*x^2 + " << par1 << "*x + " << par2 << ")" << endl;
  cout << "error par0/1/2=" << par0_err << "/" << par1_err << "/" << par2_err << endl;
  int num=1;
  //TCanvas *c2 = new TCanvas("c2","",800,600);
  /*for(int a=0;a<2;a++){
    for(int b=0;b<2;b++){
      for(int c=0;c<2;c++){
        TF1 *fix_func = new TF1(Form("fix_func%d",num),"exp([0]*x*x+[1]*x+[2])",0,fit_max);
        float fix_par0 = par0 + pm_par0_err[a];
        float fix_par1 = par1 + pm_par1_err[b];
        float fix_par2 = par2 + pm_par2_err[c];
        cout << "fix par0/1/2=" << fix_par0 << "/" << fix_par1 << "/" << fix_par2 << endl;
        fix_func->FixParameter(0,fix_par0);
        fix_func->FixParameter(1,fix_par1);
        fix_func->FixParameter(2,fix_par2);
        double total_events = fix_func->Integral(0,125)/25;
        cout << "events in SR is " << total_events << endl;
        fix_func->SetLineColor(num);
        //if(num==1) fix_func->Draw();
        fix_func->Draw("same");
        num++;
      }
    }
  }*/
  
  TF1 *fix_func = new TF1(Form("fix_func%d",num),"exp([0]*x*x+[1]*x+[2])",0,fit_max);
  fix_func->FixParameter(0,par0);
  fix_func->FixParameter(1,par1);
  fix_func->FixParameter(2,par2);
  fix_func->SetLineColor(4);
  fix_func->Draw("same");
  //c->SaveAs(Form("hist/total_distance_ntag_cut%d_nring%d_mulike%d_michel%d_fit.pdf",cut,nring,mulike,michel));
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
}

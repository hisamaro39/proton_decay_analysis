#include <vector>
void compare_sig_bkg(){
  string type[] = {"p_epi","p_mupi","p_eee","p_mumumu","p_emumu","p_mumue","p_muee","p_eemu","fcmc","fcdt"};
  int mode_id = 2;

  int nring=1;
  int mulike=0;
  int michel=0;
  string sk_period = "sk4";
  float proton_life_time = 7.9;//*10^32year

  gROOT->SetStyle("Plain");
  gStyle->SetOptStat(0);
  vector<int> dology,dorebin;
  dology.clear();dorebin.clear();
  vector<string> sig_name,sig_name_2,bkg_name,bkg_name_2,hist_name;
  sig_name.clear();sig_name_2.clear();bkg_name.clear();bkg_name_2.clear();hist_name.clear();
  
  //hist_name.push_back("nRing_cut1");
  //dology.push_back(0);dorebin.push_back(0);combine.push_back(0);

  //hist_name.push_back("nElikeRing_angle_nring3_cut1");
  //dology.push_back(0);dorebin.push_back(0);combine.push_back(0);

  //hist_name.push_back("nElikeRing_angle_nring2_cut1");
  //dology.push_back(0);dorebin.push_back(0);combine.push_back(0);
  
  //hist_name.push_back("nMulikeRing_angle_nring3_cut1");
  //dology.push_back(0);dorebin.push_back(0);combine.push_back(0);

  //hist_name.push_back("nMulikeRing_angle_nring2_cut1");
  //dology.push_back(0);dorebin.push_back(0);combine.push_back(0);

  //hist_name.push_back("n_michel_electron_cut3");
  //dology.push_back(0);dorebin.push_back(0);combine.push_back(0);

  hist_name.push_back("mom_proton_reco_masscut_cut5");
  dology.push_back(0);dorebin.push_back(5);

  hist_name.push_back("mass_proton_reco_momcut_cut5");
  dology.push_back(0);dorebin.push_back(5);

  //hist_name.push_back("ntag_multiplicity_cut4");
  //dology.push_back(0);dorebin.push_back(0);combine.push_back(0);
  
  //hist_name.push_back("mass_two_elike_reco_cut4");
  //dology.push_back(0);dorebin.push_back(5);combine.push_back(0);

  TH1 *first_hist;
  TFile *input_bkg = TFile::Open(Form("../output/fcmc.%s.mode_%s.root",sk_period.c_str(),type[mode_id].c_str()));//bkg
  TFile *input_sig = TFile::Open(Form("../output/%s_miura.%s.mode_%s_miura.root",type[mode_id].c_str(),sk_period.c_str(),type[mode_id].c_str()));//sig

  //calculate signal scale
  TH1* h_sig_cutflow = (TH1*) input_sig->Get(Form("cut_flow_nring%d_mulike%d_michel%d",nring,mulike,michel));
  float signal_scale = 0.;
  float initial_events = h_sig_cutflow->GetBinContent(1); 
  float exp_n_pdk = 640./proton_life_time;//for sk4 
  cout << "initial events is " << initial_events << endl;
  cout << "proton life time is " << proton_life_time << "*10^32year" << endl;
  cout << "expected # of proton decay is " << exp_n_pdk << endl;
  float signal_scale = exp_n_pdk/initial_events;
  cout << "signal scale is " << signal_scale << endl;

  for(int s=0;s<hist_name.size();s++){
    TH1* h_sig_all = (TH1*) input_sig->Get(Form("%s_nring%d_mulike%d_michel%d",hist_name[s].c_str(),nring,mulike,michel));
    TH1* h_sig_free = (TH1*) input_sig->Get(Form("%s_nring%d_mulike%d_michel%d_fp1",hist_name[s].c_str(),nring,mulike,michel));
    TH1* h_bkg = (TH1*) input_bkg->Get(Form("%s_nring%d_mulike%d_michel%d",hist_name[s].c_str(),nring,mulike,michel));
    h_sig_all->Rebin(dorebin[s]);
    h_sig_free->Rebin(dorebin[s]);
    h_bkg->Rebin(dorebin[s]);
    
    TCanvas *c = new TCanvas(Form("canvas%d",s),"",800,600);
    if(dology[s]) c->SetLogy();
      h_sig_all->Scale(signal_scale);
      h_sig_free->Scale(signal_scale);
      h_bkg->SetLineWidth(2);
      h_sig_all->SetLineColor(2);
      h_sig_all->SetLineWidth(2);
      h_sig_all->SetLineStyle(2);
      h_sig_free->SetLineColor(4);
      h_sig_free->SetLineWidth(2);
      h_sig_free->SetLineStyle(2);
      h_sig_all->Draw("hist E0");
      h_bkg->Draw("hist E0 same");
      h_sig_free->Draw("hist E0 same");

    stringstream str_mulike,str_michel,str_nring,str_lifetime;
    str_mulike << mulike;
    str_michel << michel;
    str_nring << nring;
    str_lifetime << proton_life_time;
    string save_name = "";
    save_name = "hist/compare_sig_bkg_" + hist_name[s];
    save_name += "_nring" + str_nring.str() + "_mulike" + str_mulike.str() +  "_michel" + str_michel.str(); 
    save_name += "_mode_" + type[mode_id] + "_" + sk_period + "_scale" + str_lifetime.str() + ".pdf"; 
    cout << "save_name = " << save_name << endl;
    c->SaveAs(save_name.c_str());
  }

}

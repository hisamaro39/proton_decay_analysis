#include <vector>
void compare_combined_sig_bkg_temp(){
  string type[] = {"p_epi","p_mupi","p_eee","p_mumumu","p_emumu","p_muee","fcmc","fcdt"};
  int mode_id = 2;
  int nring=0;//for no combine
  int mulike=0;
  int michel=0;
  string sk_period = "sk4";

  int range_momentum[] = {100,250,400,630,1000,2500,5000,10000,100000};
  int n_range_momentum = sizeof(range_momentum)/sizeof(int);

  gROOT->SetStyle("Plain");
  gStyle->SetOptStat(0);
  vector<int> dology,dorebin,combine;
  dology.clear();dorebin.clear();combine.clear();
  vector<string> sig_name,sig_name_2,bkg_name,bkg_name_2,hist_name;
  sig_name.clear();sig_name_2.clear();bkg_name.clear();bkg_name_2.clear();hist_name.clear();

  /*sig_name.push_back("distance_to_wall_cut1_nring0_mulike0_michel0");
  add_ratio.push_back(1);scale.push_back(0);dology.push_back(1);dorebin.push_back(4);

  sig_name.push_back("visible_energy_cut1_nring0_mulike0_michel0");
  add_ratio.push_back(1);scale.push_back(0);dology.push_back(1);dorebin.push_back(2);

  sig_name.push_back("nhit_OD_cut1_nring0_mulike0_michel0");
  add_ratio.push_back(1);scale.push_back(0);dology.push_back(1);dorebin.push_back(2);*/
  
  //hist_name.push_back("nRing_cut1");
  //dology.push_back(1);dorebin.push_back(0);combine.push_back(0);

  //hist_name.push_back("nElikeRing_nring3_cut1");
  //dology.push_back(1);dorebin.push_back(0);combine.push_back(0);

  //hist_name.push_back("nElikeRing_nring2_cut1");
  //dology.push_back(1);dorebin.push_back(0);combine.push_back(0);

  //hist_name.push_back("n_michel_electron_cut1");
  //dology.push_back(1);dorebin.push_back(0);combine.push_back(1);

  //hist_name.push_back("mom_proton_reco_cut3");
  //dology.push_back(1);dorebin.push_back(5);combine.push_back(1);

  hist_name.push_back("mass_proton_reco_cut4");
  dology.push_back(0);dorebin.push_back(5);combine.push_back(1);

  /*sig_name.push_back("nElikeRing_nring3_cut1_nring0_mulike0_michel0");
  add_ratio.push_back(1);scale.push_back(0);dology.push_back(1);dorebin.push_back(0);

  sig_name.push_back("nElikeRing_nring2_cut1_nring0_mulike0_michel0");
  add_ratio.push_back(1);scale.push_back(0);dology.push_back(1);dorebin.push_back(0);

  sig_name.push_back("nMulikeRing_nring3_cut1_nring0_mulike0_michel0");
  add_ratio.push_back(1);scale.push_back(0);dology.push_back(1);dorebin.push_back(0);

  sig_name.push_back("nMulikeRing_nring2_cut1_nring0_mulike0_michel0");
  add_ratio.push_back(1);scale.push_back(0);dology.push_back(1);dorebin.push_back(0);

  sig_name.push_back("n_michel_electron_cut3_nring0_mulike0_michel0");
  add_ratio.push_back(1);scale.push_back(0);dology.push_back(1);dorebin.push_back(0);

  sig_name.push_back("mass_pi0_reco_elike2_cut4_nring0_mulike0_michel0");
  add_ratio.push_back(1);scale.push_back(0);dology.push_back(1);dorebin.push_back(2);

  sig_name.push_back("ntag_multiplicity_cut5_nring0_mulike0_michel0");
  add_ratio.push_back(1);scale.push_back(0);dology.push_back(1);dorebin.push_back(0);

  sig_name.push_back("mom_proton_reco_cut6_nring0_mulike0_michel0");
  add_ratio.push_back(1);scale.push_back(0);dology.push_back(1);dorebin.push_back(5);

  sig_name.push_back("mass_proton_reco_cut6_nring0_mulike0_michel0");
  add_ratio.push_back(1);scale.push_back(0);dology.push_back(1);dorebin.push_back(5);*/

  TH1 *first_hist;
  TFile *input_bkg = TFile::Open(Form("../output/fcmc.%s.mode_%s.root",sk_period.c_str(),type[mode_id].c_str()));//bkg
  TFile *input_sig = TFile::Open(Form("../output/%s.%s.mode_%s.root",type[mode_id].c_str(),sk_period.c_str(),type[mode_id].c_str()));//sig

  for(int s=0;s<hist_name.size();s++){
    TH1* sig_fp0_1 = (TH1*) input_sig->Get(Form("%s_nring0_mulike%d_michel%d_fp0",hist_name[s].c_str(),mulike,michel));
    TH1* sig_fp0_2 = (TH1*) input_sig->Get(Form("%s_nring1_mulike%d_michel%d_fp0",hist_name[s].c_str(),mulike,michel));
    ///TH1* sig_nocomb = (TH1*) input_sig->Get(Form("%s_nring%d_mulike%d_michel%d",hist_name[s].c_str(),nring,mulike,michel));
    cout << __LINE__ << endl;
    TH1* comb_sig_fp0 = sig_fp0_1->Clone();
    comb_sig_fp0->Add(sig_fp0_2);
    TH1* sig_fp1_1 = (TH1*) input_sig->Get(Form("%s_nring0_mulike%d_michel%d_fp1",hist_name[s].c_str(),mulike,michel));
    TH1* sig_fp1_2 = (TH1*) input_sig->Get(Form("%s_nring1_mulike%d_michel%d_fp1",hist_name[s].c_str(),mulike,michel));
    //TH1* sig_fp1_nocomb = (TH1*) input_sig->Get(Form("%s_nring%d_mulike%d_michel%d_fp1",hist_name[s].c_str(),nring,mulike,michel));
    cout << __LINE__ << endl;
    TH1* comb_sig_fp1 = sig_fp1_1->Clone();
    comb_sig_fp1->Add(sig_fp1_2);
    TH1* comb_sig = comb_sig_fp0->Clone();
    comb_sig->Add(comb_sig_fp1);
    TH1* bkg_1 = (TH1*) input_bkg->Get(Form("%s_nring0_mulike%d_michel%d",hist_name[s].c_str(),mulike,michel));
    TH1* bkg_2 = (TH1*) input_bkg->Get(Form("%s_nring1_mulike%d_michel%d",hist_name[s].c_str(),mulike,michel));
    //TH1* bkg_nocomb = (TH1*) input_bkg->Get(Form("%s_nring%d_mulike%d_michel%d_fp0",hist_name[s].c_str(),nring,mulike,michel));
    cout << __LINE__ << endl;
    TH1* comb_bkg = bkg_1->Clone();
    cout << __LINE__ << endl;
    comb_bkg->Add(bkg_2);
    cout << __LINE__ << endl;
    if(dorebin[s]) {
      comb_sig->Rebin(dorebin[s]);
      comb_sig_fp1->Rebin(dorebin[s]);
      comb_bkg->Rebin(dorebin[s]);
      //sig_nocomb->Rebin(dorebin[s]);
      //sig_fp1_nocomb->Rebin(dorebin[s]);
      //bkg_nocomb->Rebin(dorebin[s]);
    }
    cout << __LINE__ << endl;
    float signal_scale = 0.7 * comb_bkg->GetMaximum()/comb_sig->GetMaximum();
    cout << __LINE__ << endl;
    //float signal_scale_nocomb = 0.5 * bkg_nocomb->GetMaximum()/sig_nocomb->GetMaximum();
    
    TCanvas *c = new TCanvas(Form("canvas%d",s),"",800,600);
    if(dology[s]) c->SetLogy();
    if(combine[s]){
      comb_sig->Scale(signal_scale);
      comb_sig_fp1->Scale(signal_scale);
      comb_bkg->SetLineWidth(2);
      comb_bkg->Draw("hist E0");
      comb_sig->SetLineColor(2);
      comb_sig->SetLineWidth(2);
      comb_sig->SetLineStyle(2);
      comb_sig->Draw("hist E0 same");
      comb_sig_fp1->SetLineColor(4);
      comb_sig_fp1->SetLineWidth(2);
      comb_sig_fp1->SetLineStyle(2);
      comb_sig_fp1->Draw("hist E0 same");
    }else{
      /*sig_nocomb->Scale(signal_scale_nocomb);
      sig_fp1_nocomb->Scale(signal_scale_nocomb);
      bkg_nocomb->SetLineWidth(2);
      bkg_nocomb->Draw("hist E0");
      sig_nocomb->SetLineColor(2);
      sig_nocomb->SetLineWidth(2);
      sig_nocomb->SetLineStyle(2);
      sig_nocomb->Draw("hist E0 same");
      sig_fp1_nocomb->SetLineColor(4);
      sig_fp1_nocomb->SetLineWidth(2);
      sig_fp1_nocomb->SetLineStyle(2);
      sig_fp1_nocomb->Draw("hist E0 same");*/
    }

    stringstream str_mulike,str_michel,str_nring;
    str_mulike << mulike;
    str_michel << michel;
    str_nring << nring;
    string save_name = "";
    save_name = "hist/compare_sig_bkg_" + hist_name[s];
    if(combine[s]) save_name += "_combine_nring0_nring1_mulike" + str_mulike.str() +  "_michel" + str_michel.str(); 
    else save_name += "_nring" + str_nring.str() + "_mulike" + str_mulike.str() +  "_michel" + str_michel.str(); 
    save_name += "_mode_" + type[mode_id] + "_" + sk_period + ".pdf"; 
    //cout << "save_name = " << save_name << endl;
    //c->SaveAs(save_name.c_str());
  }

}

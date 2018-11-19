#include <vector>
#include "TGraphAsymmErrors.h"
void make_compare_plot_th1(){
  string type[] = {"p_epi","p_mupi","p_eee","p_mumumu","p_emumu","p_muee","p_eemu","fcmc","fcdt"};

  gROOT->SetStyle("Plain");
  gStyle->SetOptStat(0);
  vector<string> hist;
  vector<vector<string>> hist_set;
  vector<vector<int>> input_type_set,mode_type_set;
  vector<int> scale,dology,input_type,add_ratio,mode_type,use_validation;
  hist.clear();hist_set.clear();scale.clear();dology.clear();input_type.clear();
  input_type_set.clear();add_ratio.clear();mode_type_set.clear();mode_type.clear();use_validation.clear();

  /*hist.push_back("nRing_cut1_nring0_mulike0_michel0");
  hist.push_back("nRing_cut1_nring0_mulike0_michel0");
  input_type.push_back(5);input_type.push_back(6);
  mode_type.push_back(5);mode_type.push_back(6); 
  add_ratio.push_back(0);scale.push_back(1);dology.push_back(0);
  hist_set.push_back(hist);hist.clear();input_type_set.push_back(input_type);input_type.clear();
  mode_type_set.push_back(mode_type);mode_type.clear();

  hist.push_back("nMulikeRing_angle_nring3_cut1_nring0_mulike0_michel0");
  hist.push_back("nMulikeRing_angle_nring3_cut1_nring0_mulike0_michel0");
  input_type.push_back(5);input_type.push_back(6);
  mode_type.push_back(5);mode_type.push_back(6); 
  add_ratio.push_back(0);scale.push_back(1);dology.push_back(0);
  hist_set.push_back(hist);hist.clear();input_type_set.push_back(input_type);input_type.clear();
  mode_type_set.push_back(mode_type);mode_type.clear();

  hist.push_back("nMulikeRing_angle_nring2_cut1_nring0_mulike0_michel0");
  hist.push_back("nMulikeRing_angle_nring2_cut1_nring0_mulike0_michel0");
  input_type.push_back(5);input_type.push_back(6);
  mode_type.push_back(5);mode_type.push_back(6); 
  add_ratio.push_back(0);scale.push_back(1);dology.push_back(0);
  hist_set.push_back(hist);hist.clear();input_type_set.push_back(input_type);input_type.clear();
  mode_type_set.push_back(mode_type);mode_type.clear();

  hist.push_back("n_michel_electron_cut1_nring0_mulike0_michel0");
  hist.push_back("n_michel_electron_cut1_nring0_mulike0_michel0");
  input_type.push_back(5);input_type.push_back(6);
  mode_type.push_back(5);mode_type.push_back(6); 
  add_ratio.push_back(0);scale.push_back(1);dology.push_back(0);
  hist_set.push_back(hist);hist.clear();input_type_set.push_back(input_type);input_type.clear();
  mode_type_set.push_back(mode_type);mode_type.clear();

  hist.push_back("mass_proton_reco_cut3_nring2_mulike1_michel1");
  hist.push_back("mass_proton_reco_cut3_nring2_mulike1_michel1");
  input_type.push_back(5);input_type.push_back(6);
  mode_type.push_back(5);mode_type.push_back(6); 
  add_ratio.push_back(0);scale.push_back(1);dology.push_back(0);
  hist_set.push_back(hist);hist.clear();input_type_set.push_back(input_type);input_type.clear();
  mode_type_set.push_back(mode_type);mode_type.clear();

  hist.push_back("mom_proton_reco_cut3_nring2_mulike1_michel1");
  hist.push_back("mom_proton_reco_cut3_nring2_mulike1_michel1");
  input_type.push_back(5);input_type.push_back(6);
  mode_type.push_back(5);mode_type.push_back(6); 
  add_ratio.push_back(0);scale.push_back(1);dology.push_back(0);
  hist_set.push_back(hist);hist.clear();input_type_set.push_back(input_type);input_type.clear();
  mode_type_set.push_back(mode_type);mode_type.clear();

  hist.push_back("true_mom_muon");
  hist.push_back("true_mom_electron");
  hist.push_back("true_mom_muon");
  hist.push_back("true_mom_electron");
  input_type.push_back(5);input_type.push_back(5);input_type.push_back(6);input_type.push_back(6);
  mode_type.push_back(5);mode_type.push_back(5);mode_type.push_back(6);mode_type.push_back(6);
  add_ratio.push_back(0);scale.push_back(0);dology.push_back(0);use_validation.push_back(1);
  hist_set.push_back(hist);hist.clear();input_type_set.push_back(input_type);input_type.clear();
  mode_type_set.push_back(mode_type);mode_type.clear();

  hist.push_back("n_true_decayE");
  hist.push_back("n_true_decayE");
  input_type.push_back(5);input_type.push_back(6);  
  mode_type.push_back(5);mode_type.push_back(6);
  add_ratio.push_back(0);scale.push_back(0);dology.push_back(0);use_validation.push_back(1);
  hist_set.push_back(hist);hist.clear();input_type_set.push_back(input_type);input_type.clear();
  mode_type_set.push_back(mode_type);mode_type.clear();

  hist.push_back("nMulikeRing_angle_nring3_trueDecayE0");
  hist.push_back("nMulikeRing_angle_nring3_trueDecayE1");
  input_type.push_back(6);input_type.push_back(6);  
  mode_type.push_back(6);mode_type.push_back(6);
  add_ratio.push_back(0);scale.push_back(1);dology.push_back(0);use_validation.push_back(1);
  hist_set.push_back(hist);hist.clear();input_type_set.push_back(input_type);input_type.clear();
  mode_type_set.push_back(mode_type);mode_type.clear();*/

  /*hist.push_back("true_mom_1st_lepton_fp0");
  hist.push_back("true_mom_2nd_lepton_fp0");
  hist.push_back("true_mom_3rd_lepton_fp0");
  input_type.push_back(5);input_type.push_back(5);input_type.push_back(5);  
  mode_type.push_back(5);mode_type.push_back(5);mode_type.push_back(5);
  add_ratio.push_back(0);scale.push_back(0);dology.push_back(0);use_validation.push_back(1);
  hist_set.push_back(hist);hist.clear();input_type_set.push_back(input_type);input_type.clear();
  mode_type_set.push_back(mode_type);mode_type.clear();

  hist.push_back("true_mass_1st_2nd_fp1");
  hist.push_back("true_mass_2nd_3rd_fp1");
  hist.push_back("true_mass_3rd_1st_fp1");
  input_type.push_back(5);input_type.push_back(5);input_type.push_back(5);  
  mode_type.push_back(5);mode_type.push_back(5);mode_type.push_back(5);
  add_ratio.push_back(0);scale.push_back(0);dology.push_back(0);use_validation.push_back(1);
  hist_set.push_back(hist);hist.clear();input_type_set.push_back(input_type);input_type.clear();
  mode_type_set.push_back(mode_type);mode_type.clear();*/

  /*hist.push_back("residual_total_mass_nring3_mulike0");
  hist.push_back("residual_total_mass_nring3_mulike1");
  input_type.push_back(2);input_type.push_back(5);
  mode_type.push_back(2);mode_type.push_back(5);
  add_ratio.push_back(0);scale.push_back(1);dology.push_back(0);
  hist_set.push_back(hist);hist.clear();input_type_set.push_back(input_type);input_type.clear();
  mode_type_set.push_back(mode_type);mode_type.clear();

  hist.push_back("residual_total_gamma_mass_nring3_mulike0");
  hist.push_back("residual_total_mass_nring3_mulike0");
  input_type.push_back(0);input_type.push_back(2);
  mode_type.push_back(0);mode_type.push_back(2);
  add_ratio.push_back(0);scale.push_back(1);dology.push_back(0);
  hist_set.push_back(hist);hist.clear();input_type_set.push_back(input_type);input_type.clear();
  mode_type_set.push_back(mode_type);mode_type.clear();

  hist.push_back("residual_total_gamma_mass_nring3_mulike1");
  hist.push_back("residual_total_mass_nring3_mulike1");
  input_type.push_back(1);input_type.push_back(5);
  mode_type.push_back(1);mode_type.push_back(5);
  add_ratio.push_back(0);scale.push_back(1);dology.push_back(0);
  hist_set.push_back(hist);hist.clear();input_type_set.push_back(input_type);input_type.clear();
  mode_type_set.push_back(mode_type);mode_type.clear();

  hist.push_back("diff_opening_angle_muon_mom200_300");
  hist.push_back("diff_opening_angle_muon_mom300_400");
  hist.push_back("diff_opening_angle_muon_mom400_500");
  input_type.push_back(5);input_type.push_back(5);input_type.push_back(5);
  mode_type.push_back(5);mode_type.push_back(5);mode_type.push_back(5);
  add_ratio.push_back(0);scale.push_back(1);dology.push_back(0);
  hist_set.push_back(hist);hist.clear();input_type_set.push_back(input_type);input_type.clear();
  mode_type_set.push_back(mode_type);mode_type.clear();

  hist.push_back("true_mom_muon");
  hist.push_back("true_mom_gamma");
  input_type.push_back(1);input_type.push_back(1);
  mode_type.push_back(1);mode_type.push_back(1);
  add_ratio.push_back(0);scale.push_back(0);dology.push_back(0);
  hist_set.push_back(hist);hist.clear();input_type_set.push_back(input_type);input_type.clear();
  mode_type_set.push_back(mode_type);mode_type.clear();

  hist.push_back("true_min_mom_lepton");
  hist.push_back("true_mid_mom_lepton");
  hist.push_back("true_max_mom_lepton");
  input_type.push_back(3);input_type.push_back(3);input_type.push_back(3);
  mode_type.push_back(3);mode_type.push_back(3); mode_type.push_back(3);
  add_ratio.push_back(0);scale.push_back(0);dology.push_back(0);
  hist_set.push_back(hist);hist.clear();input_type_set.push_back(input_type);input_type.clear();
  mode_type_set.push_back(mode_type);mode_type.clear();

  hist.push_back("true_min_mom_lepton_nring2");
  hist.push_back("true_mid_mom_lepton_nring2");
  hist.push_back("true_max_mom_lepton_nring2");
  input_type.push_back(3);input_type.push_back(3);input_type.push_back(3);
  mode_type.push_back(3);mode_type.push_back(3); mode_type.push_back(3);
  add_ratio.push_back(0);scale.push_back(0);dology.push_back(0);
  hist_set.push_back(hist);hist.clear();input_type_set.push_back(input_type);input_type.clear();
  mode_type_set.push_back(mode_type);mode_type.clear();*/

  /*hist.push_back("true_min_mom_lepton_nring4");
  hist.push_back("true_mid_mom_lepton_nring4");
  hist.push_back("true_max_mom_lepton_nring4");
  input_type.push_back(3);input_type.push_back(3);input_type.push_back(3);
  mode_type.push_back(3);mode_type.push_back(3); mode_type.push_back(3);
  add_ratio.push_back(0);scale.push_back(0);dology.push_back(0);
  hist_set.push_back(hist);hist.clear();input_type_set.push_back(input_type);input_type.clear();
  mode_type_set.push_back(mode_type);mode_type.clear();

  hist.push_back("true_min_angle_lepton_lepton_nring2");
  hist.push_back("true_min_angle_lepton_lepton_nring3");
  input_type.push_back(3);input_type.push_back(3);
  mode_type.push_back(3);mode_type.push_back(3);
  add_ratio.push_back(0);scale.push_back(0);dology.push_back(0);
  hist_set.push_back(hist);hist.clear();input_type_set.push_back(input_type);input_type.clear();
  mode_type_set.push_back(mode_type);mode_type.clear();

  hist.push_back("nDecayE_true_decayE_1");
  hist.push_back("nDecayE_true_decayE_2");
  input_type.push_back(4);input_type.push_back(4);
  mode_type.push_back(4);mode_type.push_back(4);
  add_ratio.push_back(0);scale.push_back(0);dology.push_back(0);
  hist_set.push_back(hist);hist.clear();input_type_set.push_back(input_type);input_type.clear();
  mode_type_set.push_back(mode_type);mode_type.clear();*/

  /*hist.push_back("true_angle_min_mid_lepton_nring4");
  hist.push_back("true_angle_min_max_lepton_nring4");
  hist.push_back("true_angle_mid_max_lepton_nring4");
  input_type.push_back(2);input_type.push_back(2);input_type.push_back(2);
  mode_type.push_back(2);mode_type.push_back(2); mode_type.push_back(2);
  add_ratio.push_back(0);scale.push_back(0);dology.push_back(0);
  hist_set.push_back(hist);hist.clear();input_type_set.push_back(input_type);input_type.clear();
  mode_type_set.push_back(mode_type);mode_type.clear();

  hist.push_back("true_min_angle_lepton_lepton_nring2");
  hist.push_back("true_min_angle_lepton_lepton_nring3");
  input_type.push_back(2);input_type.push_back(2);
  mode_type.push_back(2);mode_type.push_back(2);
  add_ratio.push_back(0);scale.push_back(0);dology.push_back(0);
  hist_set.push_back(hist);hist.clear();input_type_set.push_back(input_type);input_type.clear();
  mode_type_set.push_back(mode_type);mode_type.clear();

  hist.push_back("nRing_cut1_nring0_mulike0_michel0_fp0");
  hist.push_back("nRing_cut1_nring0_mulike0_michel0");
  hist.push_back("nRing_cut1_nring0_mulike0_michel0_fp1");
  input_type.push_back(6);input_type.push_back(3);input_type.push_back(3);
  mode_type.push_back(3);mode_type.push_back(3); mode_type.push_back(3);
  add_ratio.push_back(0);scale.push_back(0);dology.push_back(1);
  hist_set.push_back(hist);hist.clear();input_type_set.push_back(input_type);input_type.clear();
  mode_type_set.push_back(mode_type);mode_type.clear();*/

  /*hist.push_back("diff_opening_angle_diff_vertex_r0_5");
  hist.push_back("diff_opening_angle_diff_vertex_r5_10");
  hist.push_back("diff_opening_angle_diff_vertex_r10_15");
  hist.push_back("diff_opening_angle_diff_vertex_r15_20");
  hist.push_back("diff_opening_angle_diff_vertex_r20_25");
  hist.push_back("diff_opening_angle_diff_vertex_r25_30");
  input_type.push_back(3);input_type.push_back(3);input_type.push_back(3);
  input_type.push_back(3);input_type.push_back(3);input_type.push_back(3);
  mode_type.push_back(3);mode_type.push_back(3); mode_type.push_back(3);
  mode_type.push_back(3);mode_type.push_back(3); mode_type.push_back(3);
  add_ratio.push_back(0);scale.push_back(1);dology.push_back(1);
  hist_set.push_back(hist);hist.clear();input_type_set.push_back(input_type);input_type.clear();
  mode_type_set.push_back(mode_type);mode_type.clear();

  hist.push_back("diff_opening_angle_mom200_300");
  hist.push_back("diff_opening_angle_mom300_400");
  hist.push_back("diff_opening_angle_mom400_500");
  input_type.push_back(3);input_type.push_back(3);input_type.push_back(3);
  mode_type.push_back(3);mode_type.push_back(3); mode_type.push_back(3);
  add_ratio.push_back(0);scale.push_back(1);dology.push_back(1);
  hist_set.push_back(hist);hist.clear();input_type_set.push_back(input_type);input_type.clear();
  mode_type_set.push_back(mode_type);mode_type.clear();*/

  /*hist.push_back("residual_mmom_mom200_300");
  hist.push_back("residual_mmom_mom300_400");
  hist.push_back("residual_mmom_mom400_500");
  input_type.push_back(3);input_type.push_back(3);input_type.push_back(3);
  mode_type.push_back(3);mode_type.push_back(3); mode_type.push_back(3);
  add_ratio.push_back(0);scale.push_back(1);dology.push_back(0);
  hist_set.push_back(hist);hist.clear();input_type_set.push_back(input_type);input_type.clear();
  mode_type_set.push_back(mode_type);mode_type.clear();*/

  /*hist.push_back("residual_opening_angle_angle_mulike");
  hist.push_back("residual_opening_angle_angle_elike");
  input_type.push_back(3);input_type.push_back(3);
  mode_type.push_back(3);mode_type.push_back(3);;
  add_ratio.push_back(0);scale.push_back(1);dology.push_back(0);
  hist_set.push_back(hist);hist.clear();input_type_set.push_back(input_type);input_type.clear();
  mode_type_set.push_back(mode_type);mode_type.clear();

  hist.push_back("residual_opening_angle_charge_mulike");
  hist.push_back("residual_opening_angle_charge_elike");
  input_type.push_back(3);input_type.push_back(3);
  mode_type.push_back(3);mode_type.push_back(3);;
  add_ratio.push_back(0);scale.push_back(1);dology.push_back(0);
  hist_set.push_back(hist);hist.clear();input_type_set.push_back(input_type);input_type.clear();
  mode_type_set.push_back(mode_type);mode_type.clear();*/

  /*hist.push_back("n_michel_electron_cut3_nring1_mulike3_michel3_fp0");
  hist.push_back("n_michel_electron_cut3_nring1_mulike3_michel3");
  input_type.push_back(6);input_type.push_back(3);
  mode_type.push_back(3);mode_type.push_back(3); 
  add_ratio.push_back(0);scale.push_back(0);
  hist_set.push_back(hist);dology.push_back(1);hist.clear();input_type_set.push_back(input_type);input_type.clear();
  mode_type_set.push_back(mode_type);mode_type.clear();

  hist.push_back("mom_all_ring_reco_nring2_cut3_nring0_mulike0_michel0");
  hist.push_back("mom_all_ring_reco_nring2_cut3_nring0_mulike1_michel0");
  input_type.push_back(5);input_type.push_back(5);
  mode_type.push_back(5);mode_type.push_back(5);
  add_ratio.push_back(0);scale.push_back(0);
  hist_set.push_back(hist);dology.push_back(1);hist.clear();input_type_set.push_back(input_type);input_type.clear();
  mode_type_set.push_back(mode_type);mode_type.clear();

  hist.push_back("nRing_weight_osc_cut1_nring0_mulike0_michel0");
  hist.push_back("nRing_weight_osc_cut1_nring0_mulike0_michel0");
  input_type.push_back(6);input_type.push_back(5);
  mode_type.push_back(5);mode_type.push_back(5);
  add_ratio.push_back(0);scale.push_back(0);
  hist_set.push_back(hist);dology.push_back(1);hist.clear();input_type_set.push_back(input_type);input_type.clear();
  mode_type_set.push_back(mode_type);mode_type.clear();

  hist.push_back("nMulikeRing_nring3_weight_osc_cut1_nring0_mulike0_michel0");
  hist.push_back("nMulikeRing_nring3_weight_osc_cut1_nring0_mulike0_michel0");
  input_type.push_back(6);input_type.push_back(5);
  mode_type.push_back(5);mode_type.push_back(5);
  add_ratio.push_back(0);scale.push_back(0);
  hist_set.push_back(hist);dology.push_back(1);hist.clear();input_type_set.push_back(input_type);input_type.clear();
  mode_type_set.push_back(mode_type);mode_type.clear();

  hist.push_back("nMulikeRing_nring2_weight_osc_cut1_nring0_mulike0_michel0");
  hist.push_back("nMulikeRing_nring2_weight_osc_cut1_nring0_mulike0_michel0");
  input_type.push_back(6);input_type.push_back(5);
  mode_type.push_back(5);mode_type.push_back(5);
  add_ratio.push_back(0);scale.push_back(0);
  hist_set.push_back(hist);dology.push_back(1);hist.clear();input_type_set.push_back(input_type);input_type.clear();
  mode_type_set.push_back(mode_type);mode_type.clear();

  hist.push_back("n_michel_electron_weight_osc_cut3_nring1_mulike1_michel0");
  hist.push_back("n_michel_electron_weight_osc_cut3_nring1_mulike1_michel0");
  input_type.push_back(6);input_type.push_back(5);
  mode_type.push_back(5);mode_type.push_back(5);
  add_ratio.push_back(0);scale.push_back(0);
  hist_set.push_back(hist);dology.push_back(1);hist.clear();input_type_set.push_back(input_type);input_type.clear();
  mode_type_set.push_back(mode_type);mode_type.clear();

  hist.push_back("n_michel_electron_weight_osc_cut3_nring0_mulike1_michel0");
  hist.push_back("n_michel_electron_weight_osc_cut3_nring0_mulike1_michel0");
  input_type.push_back(6);input_type.push_back(5);
  mode_type.push_back(5);mode_type.push_back(5);
  add_ratio.push_back(0);scale.push_back(0);
  hist_set.push_back(hist);dology.push_back(1);hist.clear();input_type_set.push_back(input_type);input_type.clear();
  mode_type_set.push_back(mode_type);mode_type.clear();*/

  hist.push_back("mass_proton_reco_cut3_nring2_mulike0_michel0_type0_mode_pos22");
  hist.push_back("mass_proton_reco_cut3_nring2_mulike1_michel1_type1_mode_pos22");
  input_type.push_back(7);input_type.push_back(7);
  mode_type.push_back(2);mode_type.push_back(5);
  add_ratio.push_back(0);scale.push_back(0);use_validation.push_back(0);
  hist_set.push_back(hist);dology.push_back(0);hist.clear();input_type_set.push_back(input_type);input_type.clear();
  mode_type_set.push_back(mode_type);mode_type.clear();

  /*hist.push_back("mass_proton_reco_cut4_nring2_mulike0_michel0");
  hist.push_back("mass_proton_reco_cut4_nring2_mulike1_michel1");
  input_type.push_back(7);input_type.push_back(7);
  mode_type.push_back(2);mode_type.push_back(5);
  add_ratio.push_back(0);scale.push_back(0);use_validation.push_back(0);
  hist_set.push_back(hist);dology.push_back(1);hist.clear();input_type_set.push_back(input_type);input_type.clear();
  mode_type_set.push_back(mode_type);mode_type.clear();

  hist.push_back("mom_proton_reco_cut4_nring2_mulike0_michel0");
  hist.push_back("mom_proton_reco_cut4_nring2_mulike1_michel1");
  input_type.push_back(7);input_type.push_back(7);
  mode_type.push_back(2);mode_type.push_back(5);
  add_ratio.push_back(0);scale.push_back(0);use_validation.push_back(0);
  hist_set.push_back(hist);dology.push_back(1);hist.clear();input_type_set.push_back(input_type);input_type.clear();
  mode_type_set.push_back(mode_type);mode_type.clear();

  hist.push_back("mom_proton_reco_cut3_nring1_mulike3_michel0");
  hist.push_back("mom_proton_reco_weight_osc_cut3_nring1_mulike3_michel0");
  input_type.push_back(6);input_type.push_back(3);
  mode_type.push_back(3);mode_type.push_back(3);
  add_ratio.push_back(0);scale.push_back(0);
  hist_set.push_back(hist);dology.push_back(1);hist.clear();input_type_set.push_back(input_type);input_type.clear();
  mode_type_set.push_back(mode_type);mode_type.clear();

  hist.push_back("mom_proton_reco_cut5_nring0_mulike0_michel0");
  hist.push_back("mom_proton_reco_cut5_nring1_mulike0_michel0");
  hist.push_back("mom_proton_reco_weight_osc_cut5_nring0_mulike0_michel0");
  hist.push_back("mom_proton_reco_weight_osc_cut5_nring1_mulike0_michel0");
  input_type.push_back(2);input_type.push_back(2);input_type.push_back(6);input_type.push_back(6);
  mode_type.push_back(2);mode_type.push_back(2);mode_type.push_back(2);mode_type.push_back(2);
  add_ratio.push_back(0);scale.push_back(0);
  hist_set.push_back(hist);dology.push_back(1);hist.clear();input_type_set.push_back(input_type);input_type.clear();
  mode_type_set.push_back(mode_type);mode_type.clear();*/

  TFile *input;
  TH1 *first_hist;
  for(int s=0;s<hist_set.size();s++){
    float evtmax=0,evtmax_scale=0,ratio_max=0,ratio_min=99999,evt_min,evtmin_scale=1;
    float very_small_evtmax = 99;
    for(int ss=0;ss<hist_set[s].size();ss++){//decide max event of the hist
      if(use_validation[s]) input = TFile::Open(Form("../output/%s.sk4.mode_%s_validation.root",type[input_type_set[s].at(ss)].c_str(),type[mode_type_set[s].at(ss)].c_str()));
      else input = TFile::Open(Form("../output/%s.sk4.mode_%s_check_bkg.root",type[input_type_set[s].at(ss)].c_str(),type[mode_type_set[s].at(ss)].c_str()));
      TH1* temp_hist = (TH1*) input->Get(hist_set[s].at(ss).c_str());
      if(temp_hist->GetMaximum() > evtmax) evtmax = temp_hist->GetMaximum();
      if(temp_hist->GetEntries() && temp_hist->GetMaximum()/temp_hist->Integral() > evtmax_scale){ 
        evtmax_scale = temp_hist->GetMaximum()/temp_hist->Integral();
        evtmin_scale = 1./temp_hist->Integral();
      }
      if(evtmax<1) very_small_evtmax = evtmax;
      if(ss==0) first_hist=temp_hist;
      else{
        TH1 *temp_ratio = (TH1*) temp_hist->Clone("temp_ratio");
        temp_ratio->Divide(first_hist);
        for(int n=0;n<temp_ratio->GetNbinsX();n++){
          float this_ratio = temp_ratio->GetBinContent(n+1);
          if(this_ratio<ratio_min) ratio_min=this_ratio;
          if(this_ratio>ratio_max) ratio_max=this_ratio;
        }
      }
    }
    TCanvas *c = new TCanvas(Form("canvas%d",s),"",800,600);
    if(add_ratio[s]){
      TPad* p1 = new TPad("main","main",0.0,0.2,1.0,1.0,10,0,0);
      p1->SetNumber(1);
      p1->SetBottomMargin(0);
      p1->Draw();
      TPad* p2 = new TPad("ratio","ratio",0.0,0.0,1.0,0.2,10,0,0);
      p2->SetTopMargin(0);
      p2->SetBottomMargin(0.5);
      p2->SetNumber(2);
      p2->Draw();
    }
    if(dology[s]) c->SetLogy();
    string save_name = "";
    save_name = "hist/compare_";
    for(int h=0;h<hist_set[s].size();h++){
      if(use_validation[s]) input = TFile::Open(Form("../output/%s.sk4.mode_%s_validation.root",type[input_type_set[s].at(h)].c_str(),type[mode_type_set[s].at(h)].c_str()));
      else input = TFile::Open(Form("../output/%s.sk4.mode_%s_check_bkg.root",type[input_type_set[s].at(h)].c_str(),type[mode_type_set[s].at(h)].c_str()));
      TH1* this_hist = (TH1*) input->Get(hist_set[s].at(h).c_str());
      if(scale[s] && this_hist->GetEntries()) this_hist->Scale(1./this_hist->Integral());
      this_hist->SetLineWidth(2);
      if(add_ratio[s]) c->cd(1);
      if(h==0){
        if(scale[s]) this_hist->GetYaxis()->SetRangeUser(0,evtmax_scale*1.2);
        else this_hist->GetYaxis()->SetRangeUser(0,evtmax*1.2);
        if(dology[s]){
          if(scale[s]) this_hist->GetYaxis()->SetRangeUser(evtmin_scale,evtmax_scale*1.2);
          else {
            this_hist->GetYaxis()->SetRangeUser(0.001,evtmax*1.2);
            if(very_small_evtmax<1) this_hist->GetYaxis()->SetRangeUser(very_small_evtmax*0.01,evtmax*1.2);
          }
        }
        this_hist->SetLineColor(1);
        this_hist->Draw("hist E0");
        //this_hist->Draw();
        first_hist = this_hist;
      }
      else{
        this_hist->SetLineColor(h+1);
        this_hist->Draw("same hist E0");
        TH1 *ratio_hist = (TH1*) this_hist->Clone("ratio_plot");
        ratio_hist->Divide(first_hist);
        float xmin = ratio_hist->GetBinLowEdge(1);
        float xmax = ratio_hist->GetBinLowEdge(ratio_hist->GetNbinsX())+ratio_hist->GetBinWidth(ratio_hist->GetNbinsX());
        if(add_ratio[s]){
          c->cd(2);
          if(h==1){
            TH1* frame;
            frame=gPad->DrawFrame(xmin, ratio_min*0.9, xmax, ratio_max*1.1);
            if(ratio_max<1) frame=gPad->DrawFrame(xmin, ratio_min*0.9, xmax, 1.1);
            if(ratio_min>1) frame=gPad->DrawFrame(xmin, 0.9, xmax, ratio_max*1.1);
            frame->GetYaxis()->SetLabelSize(0.1);
            frame->GetXaxis()->SetLabelSize(0.2);
          }
          ratio_hist->Draw("same");
          TLine *line = new TLine(xmin,1,xmax,1);
          line->SetLineStyle(2);
          line->Draw();
        }
      }
      save_name += hist_set[s].at(h);
      save_name += "_" + type[input_type_set[s].at(h)] + "_" + type[input_type_set[s].at(h)];
    }
    save_name += ".pdf";
    //c->SaveAs(save_name.c_str());
    //c->SaveAs("hist/temp.pdf");

  }

}

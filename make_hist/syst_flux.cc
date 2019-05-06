void syst_flux(){
  //refer atmpd/src/analysis/syserror.database/Systematic.cc
  string mode_name = "p_epi";
  int period = 4;//5:sk1-4

  TChain *ch_mc = new TChain("osc_tuple");

  cout << "mode " << mode_name << " sk" << period << endl;
  if(period==5){
    for(int p=1;p<5;p++){
      ch_mc->Add(Form("../output/fcmc_final.sk%d.mode_%s_total_box_tree.root",p,mode_name.c_str()));
    }
  }else{
    ch_mc->Add(Form("../output/fcmc_final.sk%d.mode_%s_total_box_tree.root",period,mode_name.c_str()));
  }

  int nring_mc,nmulike_mc,mode,ipnu;
  float mc_weight,osc_weight,livetime_weight,total_mass_mc,total_mom_mc;
  float msdir_x_mc[5],msdir_y_mc[5],msdir_z_mc[5],pnu;
  float emom_mc[5],error[30];
  ch_mc->SetBranchAddress("nring",&nring_mc);
  ch_mc->SetBranchAddress("nmulike",&nmulike_mc);
  ch_mc->SetBranchAddress("mc_weight",&mc_weight);
  ch_mc->SetBranchAddress("pnu",&pnu);
  ch_mc->SetBranchAddress("error",&error);
  ch_mc->SetBranchAddress("ipnu",&ipnu);
  ch_mc->SetBranchAddress("osc_weight",&osc_weight);
  ch_mc->SetBranchAddress("livetime_weight",&livetime_weight);
  ch_mc->SetBranchAddress("mode",&mode);
  ch_mc->SetBranchAddress("total_mass",&total_mass_mc);
  ch_mc->SetBranchAddress("total_mom",&total_mom_mc);
  ch_mc->SetBranchAddress("msdir_x",&msdir_x_mc);
  ch_mc->SetBranchAddress("msdir_y",&msdir_y_mc);
  ch_mc->SetBranchAddress("msdir_z",&msdir_z_mc);
  ch_mc->SetBranchAddress("emom",&emom_mc);
  TFile *output = new TFile(Form("output/output_check_excess_p_eee_sk%d.root",period),"recreate");
  string error_name[] = {"abs_norm_E_lt_1GeV","abs_norm_E_gt_1GeV","nu_nubar_ratio_E_lt_1GeV",//2
  "nu_nubar_ratio_1_E_10GeV","nu_nubar_ratio_E_gt_10GeV","nuebar_nue_E_lt_1GeV",//5
  "nuebar_nue_1_E_10GeV","nuebar_nue_E_gt_10GeV","numubar_numu_E_lt_1GeV",//8
  "numubar_numu_1_E_10GeV","numubar_numu_E_gt_10GeV","up_down_ratio","horizontal_vertical_ratio",//12
  "K_pi_ratio","nu_path","axial_mass_QE_and_1pi","CCQE_xsec_ratio","CCQE_nu_nubar_ratio",//17
  "CCQE_numu_nue_ratio","single_meson_xsec","pi0_qpi_ratio","nubar_nu_1pi_ratio",//21
  "DIS_model_difference","DIS_xsec","coherent_pi_xsec","NC_CC_ratio"};
  int num_error = sizeof(error_name)/sizeof(string);
  vector<float> total_events_up,total_error_up,total_events_down,total_error_down;
  float total_events_def=0,total_error_def=0;
  for(int i=0;i<num_error;i++){//initialize
    total_events_up.push_back(0);
    total_error_up.push_back(0);
    total_events_down.push_back(0);
    total_error_down.push_back(0);
  }
  cout << "# of error is " << num_error << endl;
  cout << "# of entries is " << ch_mc->GetEntries() << endl;
  for (int e=0;e<ch_mc->GetEntries();e++){
    //if(e%100000==0) cout << "event" << e << endl;
    ch_mc->GetEntry(e);
    cout << "ipnu/pnu=" << ipnu << "/" << pnu << endl;
    cout << "weight mc/osc/livetime=" << mc_weight << "/" << osc_weight << "/" << livetime_weight << endl;
    cout << "mode=" << mode << endl;
    cout << "DIS model error is " << error[22] << endl;
    float weight = mc_weight*osc_weight*livetime_weight;
    total_events_def+=weight;
    total_error_def+=weight*weight;

    for(int k=0;k<num_error;k++){
      total_events_up[k]+=weight*(1+error[k]);
      total_error_up[k]+=pow(weight*(1+error[k]),2);
      total_events_down[k]+=weight*(1-error[k]);
      total_error_down[k]+=pow(weight*(1-error[k]),2);
    }
  }
  total_error_def = sqrt(total_error_def);
  for(int j=0;j<num_error;j++){
    total_error_up[j] = sqrt(total_error_up[j]);
    total_error_down[j] = sqrt(total_error_down[j]);
  }
  cout << "# of events in total signal box" << endl;
  cout << "default" << endl; 
  cout << "events:" << total_events_def << "+-" << total_error_def << endl;
  float total_ratio_flux=0.,total_error_flux=0.,total_ratio_xsec=0.,total_error_xsec=0.;
  for(int er=0;er<num_error;er++){
    float ratio_up = 100.*(total_events_up[er]-total_events_def)/total_events_def;
    float ratio_down = 100.*(total_events_down[er]-total_events_def)/total_events_def;
    cout << er << ": " << error_name[er] << endl;
    cout << "Up   events:" << total_events_up[er] << "+-" << total_error_up[er] 
      << " ratio:" << ratio_up << "+-" << ratio_up*total_error_up[er]/total_events_up[er] << "%" << endl;
    cout << "Down events:" << total_events_down[er] << "+-" << total_error_down[er] 
      << " ratio:" << ratio_down << "+-" << ratio_down*total_error_down[er]/total_events_down[er] << "%" << endl;
    if(er<15){
      total_ratio_flux += ratio_up*ratio_up;
      total_error_flux += pow(ratio_up*total_error_up[er]/total_events_up[er],2);
    }else {
      total_ratio_xsec += ratio_up*ratio_up;
      total_error_xsec += pow(ratio_up*total_error_up[er]/total_events_up[er],2);
    }
  }
  total_ratio_flux = sqrt(total_ratio_flux);
  total_error_flux = sqrt(total_error_flux);
  total_ratio_xsec = sqrt(total_ratio_xsec);
  total_error_xsec = sqrt(total_error_xsec);
  cout << "Flux total ratio is " << total_ratio_flux << "+-" << total_error_flux << "%" << endl;
  cout << "Xsec total ratio is " << total_ratio_xsec << "+-" << total_error_xsec << "%" << endl;

}

void syst_flux_2(){
  //refer atmpd/src/analysis/syserror.database/Systematic.cc
  string mode_name = "p_muee";
  int period = 4;//5:sk1-4

  TChain *ch_mc = new TChain("osc_tuple");

  cout << "mode " << mode_name << " sk" << period << endl;
  if(period==5){
    for(int p=1;p<5;p++){
      ch_mc->Add(Form("../output/fcmc_real.sk%d.mode_%s_total_box_tree.root",p,mode_name.c_str()));
    }
  }else{
    ch_mc->Add(Form("../output/fcmc_real.sk%d.mode_%s_total_box_tree.root",period,mode_name.c_str()));
  }

  TChain *ch_sample = new TChain("syst_info");
  ch_sample->Add("/disk01/usr3/raw/shared/fitted.error.friends/repro.16b.fitted.error.friends/sample_error_friend.sk4.root");
  char error_name[222];
  float sigma;
  vector<string> syst_name;
  vector<float> syst_sigma;
  ch_sample->SetBranchAddress("ErrorName",&error_name);
  ch_sample->SetBranchAddress("Sigma",&sigma);
  for(int er=0;er<160;er++){
    ch_sample->GetEntry(er);
    //cout << error_name << endl;
    syst_name.push_back(error_name);
    syst_sigma.push_back(sigma);
  }
  //for(int s=0;s<syst_name.size();s++) cout << syst_name[s] << endl;

  int nring_mc,nmulike_mc,mode,ipnu;
  float mc_weight,osc_weight,livetime_weight,total_mass_mc,total_mom_mc;
  float msdir_x_mc[5],msdir_y_mc[5],msdir_z_mc[5],pnu;
  float emom_mc[5],error[30],syst[160];
  ch_mc->SetBranchAddress("nring",&nring_mc);
  ch_mc->SetBranchAddress("nmulike",&nmulike_mc);
  ch_mc->SetBranchAddress("mc_weight",&mc_weight);
  ch_mc->SetBranchAddress("pnu",&pnu);
  ch_mc->SetBranchAddress("error",&error);
  ch_mc->SetBranchAddress("syst",&syst);
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

  vector<float> total_events_syst,total_error_syst;
  float total_events_def=0,total_error_def=0;
  for(int i=0;i<160;i++){//initialize
    total_events_syst.push_back(0);
    total_error_syst.push_back(0);
  }
  //cout << "# of error is " << num_error << endl;
  cout << "# of entries is " << ch_mc->GetEntries() << endl;
  for (int e=0;e<ch_mc->GetEntries();e++){
    //cout << "event" << e << endl;
    ch_mc->GetEntry(e);
    for(int er=0;er<160;er++){
      //cout << "er" << er << " weight=" << syst[er] << endl; 
    }
    //cout << "ipnu/pnu=" << ipnu << "/" << pnu << endl;
    //cout << "weight mc/osc/livetime=" << mc_weight << "/" << osc_weight << "/" << livetime_weight << endl;
    cout << "mode=" << mode << endl;
    cout << "DIS model error is " << syst[15] << endl;
    cout << "NEUT axial mass error is " << syst[64] << endl;
    //float weight = mc_weight*osc_weight*livetime_weight;//for fcmc_final
    float weight;
    if(period==5) weight = osc_weight;//for fcmc_real
    else weight = osc_weight*livetime_weight;
    total_events_def+=weight;
    total_error_def+=weight*weight;

    for(int k=0;k<160;k++){
      total_events_syst[k]+=weight*(1+syst[k]);
      total_error_syst[k]+=pow(weight*(1+syst[k]),2);
    }
  }
  total_error_def = sqrt(total_error_def);
  for(int j=0;j<160;j++) total_error_syst[j] = sqrt(total_error_syst[j]);
  cout << "# of events in total signal box" << endl;
  cout << "default" << endl; 
  cout << "events:" << total_events_def << "+-" << total_error_def << endl;
  //float total_ratio_flux=0.,total_error_flux=0.,total_ratio_xsec=0.,total_error_xsec=0.;
  for(int er=0;er<160;er++){
    float ratio_syst = 100.*(total_events_syst[er]-total_events_def)/total_events_def;
    if(ratio_syst>0.1){
      cout << er << ": " << syst_name[er] << " Sigma:" << syst_sigma[er] << endl;
      cout << "Syst events:" << total_events_syst[er] << "+-" << total_error_syst[er] 
        << " ratio:" << ratio_syst << "+-" << ratio_syst*total_error_syst[er]/total_events_syst[er] << "%" << endl;
    }
  }

}

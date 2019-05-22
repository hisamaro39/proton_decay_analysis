void syst_fsi(){
  //refer atmpd/src/analysis/syserror.database/Systematic.cc
  string mode_name = "p_mumumu";
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

  int nring_mc,nmulike_mc,mode,ipnu;
  float mc_weight,osc_weight,livetime_weight,total_mass_mc,total_mom_mc;
  float msdir_x_mc[5],msdir_y_mc[5],msdir_z_mc[5],pnu;
  float emom_mc[5],error[30],syst[160],fsi[17];
  ch_mc->SetBranchAddress("nring",&nring_mc);
  ch_mc->SetBranchAddress("nmulike",&nmulike_mc);
  ch_mc->SetBranchAddress("mc_weight",&mc_weight);
  ch_mc->SetBranchAddress("pnu",&pnu);
  ch_mc->SetBranchAddress("error",&error);
  ch_mc->SetBranchAddress("syst",&syst);
  ch_mc->SetBranchAddress("fsi",&fsi);
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

  vector<float> total_events_fsi,total_error_fsi;
  float total_events_def=0,total_error_def=0;
  for(int i=0;i<17;i++){//initialize
    total_events_fsi.push_back(0);
    total_error_fsi.push_back(0);
  }
  //cout << "# of error is " << num_error << endl;
  cout << "# of entries is " << ch_mc->GetEntries() << endl;
  for (int e=0;e<ch_mc->GetEntries();e++){
    //cout << "event" << e << endl;
    ch_mc->GetEntry(e);
    for(int f=0;f<17;f++){
      //cout << "fsi" << f << " weight=" << fsi[f] << endl; 
    }
    //cout << "ipnu/pnu=" << ipnu << "/" << pnu << endl;
    //cout << "weight mc/osc/livetime=" << mc_weight << "/" << osc_weight << "/" << livetime_weight << endl;
    //cout << "mode=" << mode << endl;
    //cout << "DIS model error is " << error[22] << endl;
    //float weight = mc_weight*osc_weight*livetime_weight;//for fcmc_final
    float weight;
    if(period==5) weight = osc_weight;//for fcmc_real
    else weight = osc_weight*livetime_weight;
    total_events_def+=weight;
    total_error_def+=weight*weight;

    for(int f=0;f<17;f++){
      total_events_fsi[f]+=weight*fsi[f];
      total_error_fsi[f]+=pow(weight*fsi[f],2);
    }
  }
  total_error_def = sqrt(total_error_def);
  for(int j=0;j<17;j++) total_error_fsi[j] = sqrt(total_error_fsi[j]);
  cout << "# of events in total signal box" << endl;
  cout << "default events:" << total_events_def << "+-" << total_error_def << endl;
  //float total_ratio_flux=0.,total_error_flux=0.,total_ratio_xsec=0.,total_error_xsec=0.;
  float mean=0,rms=0;
  for (int k=0;k<17;k++) mean += total_events_fsi[k]/17;
  for(int er=0;er<17;er++){
    float ratio_fsi = 100.*(total_events_fsi[er]-total_events_def)/total_events_def;
    cout << "fsi" << er <<  " events:" << total_events_fsi[er] << "+-" << total_error_fsi[er] 
      << " ratio:" << ratio_fsi << "+-" << ratio_fsi*total_error_fsi[er]/total_events_fsi[er] << "%" << endl;
    rms += pow((mean - total_events_fsi[er]),2);
  }
  rms = sqrt(rms)/17;
  cout << "mean=" << mean << endl;
  cout << "rms=" << rms << endl;
  cout << "rms/mean=" << rms/mean << endl;

}

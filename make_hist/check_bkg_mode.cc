void check_bkg_mode(){
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

  cout << "# of entries is " << ch_mc->GetEntries() << endl;
  float total=0,ccqe=0,other=0,cc_single_pi=0,cc_multi_pi=0,cc_eta=0;
  float nc_multi_pi=0,nc_dis=0,cc_dis=0,nc_single_pi=0,nc_eta=0,cc_k=0;
  for (int e=0;e<ch_mc->GetEntries();e++){
    ch_mc->GetEntry(e);
    float weight;
    cout << "ipnu/mode=" << ipnu << "/" << mode << endl;
    if(period==5) weight = osc_weight;//for fcmc_real
    else weight = osc_weight*livetime_weight;
    total+=weight;
    if(abs(mode)==1) ccqe+=weight;
    else if(abs(mode)>=11 && abs(mode)<=16) cc_single_pi+=weight;
    else if(abs(mode)==21) cc_multi_pi+=weight;
    else if(abs(mode)==22) cc_eta+=weight;
    else if(abs(mode)==23) cc_k+=weight;
    else if(abs(mode)==26) cc_dis+=weight;
    else if(abs(mode)>=31 && abs(mode)<=38) nc_single_pi+=weight;
    else if(abs(mode)==41) nc_multi_pi+=weight;
    else if(abs(mode)==42 || abs(mode)==43) nc_eta+=weight;
    else if(abs(mode)==46) nc_dis+=weight;
    else other+=weight;
  }
  cout << "total events is " << total << endl;
  cout << "ccqe events/ratio=" << ccqe << "/" << ccqe/total << endl;
  cout << "cc_single_pi events/ratio=" << cc_single_pi << "/" << cc_single_pi/total << endl;
  cout << "cc_multi_pi events/ratio=" << cc_multi_pi << "/" << cc_multi_pi/total << endl;
  cout << "cc_eta events/ratio=" << cc_eta << "/" << cc_eta/total << endl;
  cout << "cc_k events/ratio=" << cc_k << "/" << cc_k/total << endl;
  cout << "cc_dis events/ratio=" << cc_dis << "/" << cc_dis/total << endl;
  cout << "nc_single_pi events/ratio=" << nc_single_pi << "/" << nc_single_pi/total << endl;
  cout << "nc_multi_pi events/ratio=" << nc_multi_pi << "/" << nc_multi_pi/total << endl;
  cout << "nc_eta events/ratio=" << nc_eta << "/" << nc_eta/total << endl;
  cout << "nc_dis events/ratio=" << nc_dis << "/" << nc_dis/total << endl;
  cout << "other events/ratio=" << other << "/" << other/total << endl;

}

void check_excess(){
  int period = 5;//5:sk1-4

  TChain *ch_dt = new TChain("osc_tuple");
  TChain *ch_mc = new TChain("osc_tuple");

  if(period==5){
    for(int p=1;p<5;p++){
      ch_dt->Add(Form("../output/fcdt_final.sk%d.mode_p_eee_CR_tree.root",p));
      ch_mc->Add(Form("../output/fcmc_final.sk%d.mode_p_eee_CR_tree.root",p));
    }
  }else{
    ch_dt->Add(Form("../output/fcdt_final.sk%d.mode_p_eee_CR_tree.root",period));
    ch_mc->Add(Form("../output/fcmc_final.sk%d.mode_p_eee_CR_tree.root",period));
  }

  int nring_dt,nring_mc,nmulike_dt,nmulike_mc,mode;
  float weight,total_mass_dt,total_mass_mc,total_mom_dt,total_mom_mc;
  float msdir_x_dt[5],msdir_x_mc[5],msdir_y_dt[5],msdir_y_mc[5],msdir_z_dt[5],msdir_z_mc[5];
  float emom_dt[5],emom_mc[5];
  ch_dt->SetBranchAddress("nring",&nring_dt);
  ch_mc->SetBranchAddress("nring",&nring_mc);
  ch_dt->SetBranchAddress("nmulike",&nmulike_dt);
  ch_mc->SetBranchAddress("nmulike",&nmulike_mc);
  ch_mc->SetBranchAddress("weight",&weight);
  ch_mc->SetBranchAddress("mode",&mode);
  ch_dt->SetBranchAddress("total_mass",&total_mass_dt);
  ch_mc->SetBranchAddress("total_mass",&total_mass_mc);
  ch_dt->SetBranchAddress("total_mom",&total_mom_dt);
  ch_mc->SetBranchAddress("total_mom",&total_mom_mc);
  ch_dt->SetBranchAddress("msdir_x",&msdir_x_dt);
  ch_mc->SetBranchAddress("msdir_x",&msdir_x_mc);
  ch_dt->SetBranchAddress("msdir_y",&msdir_y_dt);
  ch_mc->SetBranchAddress("msdir_y",&msdir_y_mc);
  ch_dt->SetBranchAddress("msdir_z",&msdir_z_dt);
  ch_mc->SetBranchAddress("msdir_z",&msdir_z_mc);
  ch_dt->SetBranchAddress("emom",&emom_dt);
  ch_mc->SetBranchAddress("emom",&emom_mc);
  TFile *output = new TFile(Form("output/output_check_excess_p_eee_sk%d.root",period),"recreate");

  TH1 *h_nring_dt = new TH1F("h_nring_dt",";;",6,0,6);
  TH1 *h_nring_mc = new TH1F("h_nring_mc",";;",6,0,6);
  TH1 *h_total_mass_dt = new TH1F("h_total_mass_dt",";;",40,0,2000);
  TH1 *h_total_mass_mc = new TH1F("h_total_mass_mc",";;",40,0,2000);
  TH1 *h_total_mom_dt = new TH1F("h_total_mom_dt",";;",200,0,2000);
  TH1 *h_total_mom_mc = new TH1F("h_total_mom_mc",";;",200,0,2000);
  TH1 *h_total_mom_mass750_900_dt = new TH1F("h_total_mom_mass750_900_dt",";;",300,0,3000);
  TH1 *h_total_mom_mass750_900_mc = new TH1F("h_total_mom_mass750_900_mc",";;",300,0,3000);
  TH1 *h_mass12_dt = new TH1F("h_mass12_dt",";;",20,0,2000);
  TH1 *h_mass12_mc = new TH1F("h_mass12_mc",";;",20,0,2000);
  TH1 *h_mass12_cc_eta_mc = new TH1F("h_mass12_cc_eta_mc",";;",20,0,2000);
  TH1 *h_mass12_cc_dis_mc = new TH1F("h_mass12_cc_dis_mc",";;",20,0,2000);
  TH1 *h_mass12_nc_dis_mc = new TH1F("h_mass12_nc_dis_mc",";;",20,0,2000);
  TH1 *h_mass12_other_mc = new TH1F("h_mass12_other_mc",";;",20,0,2000);
  TH1 *h_mass12_cc_multi_pi_mc = new TH1F("h_mass12_cc_multi_pi_mc",";;",20,0,2000);
  TH1 *h_mass12_nc_multi_pi_mc = new TH1F("h_mass12_nc_multi_pi_mc",";;",20,0,2000);
  TH1 *h_mass12_cc_single_pi_mc = new TH1F("h_mass12_cc_single_pi_mc",";;",20,0,2000);
  TH1 *h_mass23_dt = new TH1F("h_mass23_dt",";;",20,0,2000);
  TH1 *h_mass23_mc = new TH1F("h_mass23_mc",";;",20,0,2000);
  TH1 *h_mass23_cc_eta_mc = new TH1F("h_mass23_cc_eta_mc",";;",20,0,2000);
  TH1 *h_mass23_cc_dis_mc = new TH1F("h_mass23_cc_dis_mc",";;",20,0,2000);
  TH1 *h_mass23_nc_dis_mc = new TH1F("h_mass23_nc_dis_mc",";;",20,0,2000);
  TH1 *h_mass23_other_mc = new TH1F("h_mass23_other_mc",";;",20,0,2000);
  TH1 *h_mass23_cc_multi_pi_mc = new TH1F("h_mass23_cc_multi_pi_mc",";;",20,0,2000);
  TH1 *h_mass23_nc_multi_pi_mc = new TH1F("h_mass23_nc_multi_pi_mc",";;",20,0,2000);
  TH1 *h_mass23_cc_single_pi_mc = new TH1F("h_mass23_cc_single_pi_mc",";;",20,0,2000);
  TH1 *h_mass13_dt = new TH1F("h_mass13_dt",";;",20,0,2000);
  TH1 *h_mass13_mc = new TH1F("h_mass13_mc",";;",20,0,2000);
  TH1 *h_mass13_cc_eta_mc = new TH1F("h_mass13_cc_eta_mc",";;",20,0,2000);
  TH1 *h_mass13_cc_dis_mc = new TH1F("h_mass13_cc_dis_mc",";;",20,0,2000);
  TH1 *h_mass13_nc_dis_mc = new TH1F("h_mass13_nc_dis_mc",";;",20,0,2000);
  TH1 *h_mass13_other_mc = new TH1F("h_mass13_other_mc",";;",20,0,2000);
  TH1 *h_mass13_cc_multi_pi_mc = new TH1F("h_mass13_cc_multi_pi_mc",";;",20,0,2000);
  TH1 *h_mass13_nc_multi_pi_mc = new TH1F("h_mass13_nc_multi_pi_mc",";;",20,0,2000);
  TH1 *h_mass13_cc_single_pi_mc = new TH1F("h_mass13_cc_single_pi_mc",";;",20,0,2000);
  TH1 *h_mode_mc = new TH1F("h_mode_mc",";;",120,-60,60);
  TH1 *h_mode_weight_mc = new TH1F("h_mode_weight_mc",";;",120,-60,60);

  cout << "Data" << endl;
  cout << "# of entries is " << ch_dt->GetEntries() << endl;
  TLorentzVector ring1,ring2,ring3;
  float evt_mc_3bin=0.,evt_dt_3bin=0.,evt_mc_box=0.,evt_dt_box=0.;
  for (int e=0;e<ch_dt->GetEntries();e++){
    //if(e%100000==0) cout << "event" << e << endl;
    ch_dt->GetEntry(e);
    float min_mom=999999,mid_mom=999999,max_mom=999999;
    for(int r=0;r<3;r++){
      float px = emom_dt[r]*msdir_x_dt[r];
      float py = emom_dt[r]*msdir_y_dt[r];
      float pz = emom_dt[r]*msdir_z_dt[r];
      float energy = sqrt(emom_dt[r]*emom_dt[r]+0.511*0.511);
      if(emom_dt[r]<min_mom) {
        max_mom=mid_mom;
        mid_mom=min_mom;
        min_mom=emom_dt[r];
        ring3 = ring2;
        ring2 = ring1;
        ring1.SetPxPyPzE(px,py,pz,energy);
      }
      else if(emom_dt[r]<mid_mom) {
        max_mom=mid_mom;
        mid_mom=emom_dt[r];
        ring3 = ring2;
        ring2.SetPxPyPzE(px,py,pz,energy);
      }
      else if(emom_dt[r]<max_mom) {
        max_mom=emom_dt[r];
        ring3.SetPxPyPzE(px,py,pz,energy);
      }
    }
    float mass12 = (ring1+ring2).M();
    float mass23 = (ring2+ring3).M();
    float mass13 = (ring1+ring3).M();
    h_nring_dt->Fill(nring_dt);
    h_total_mass_dt->Fill(total_mass_dt);
    h_total_mom_dt->Fill(total_mom_dt);
    if(total_mass_dt>750 && total_mass_dt<900){
      h_total_mom_mass750_900_dt->Fill(total_mom_dt);
      h_mass12_dt->Fill(mass12);
      h_mass23_dt->Fill(mass23);
      h_mass13_dt->Fill(mass13);
      evt_dt_3bin++;
    }
    if(total_mass_dt>800 && total_mass_dt<1050) evt_dt_box++;
  }

  cout << "MC" << endl;
  cout << "# of entries is " << ch_mc->GetEntries() << endl;
  for (int e=0;e<ch_mc->GetEntries();e++){
    //if(e%100000==0) cout << "event" << e << endl;
    ch_mc->GetEntry(e);
    float min_mom=999999,mid_mom=999999,max_mom=999999;
    for(int r=0;r<3;r++){
      float px = emom_mc[r]*msdir_x_mc[r];
      float py = emom_mc[r]*msdir_y_mc[r];
      float pz = emom_mc[r]*msdir_z_mc[r];
      float energy = sqrt(emom_mc[r]*emom_mc[r]+0.511*0.511);
      if(emom_mc[r]<min_mom) {
        max_mom=mid_mom;
        mid_mom=min_mom;
        min_mom=emom_mc[r];
        ring3 = ring2;
        ring2 = ring1;
        ring1.SetPxPyPzE(px,py,pz,energy);
      }
      else if(emom_mc[r]<mid_mom) {
        max_mom=mid_mom;
        mid_mom=emom_mc[r];
        ring3 = ring2;
        ring2.SetPxPyPzE(px,py,pz,energy);
      }
      else if(emom_mc[r]<max_mom) {
        max_mom=emom_mc[r];
        ring3.SetPxPyPzE(px,py,pz,energy);
      }
    }
    float mass12 = (ring1+ring2).M();
    float mass23 = (ring2+ring3).M();
    float mass13 = (ring1+ring3).M();
    h_nring_mc->Fill(nring_mc,weight);
    h_total_mass_mc->Fill(total_mass_mc,weight);
    h_total_mom_mc->Fill(total_mom_mc,weight);
    if(total_mass_mc>750 && total_mass_mc<900){
      h_total_mom_mass750_900_mc->Fill(total_mom_mc,weight);
      h_mode_mc->Fill(mode);
      h_mode_weight_mc->Fill(mode,weight);
      h_mass12_mc->Fill(mass12,weight);
      h_mass23_mc->Fill(mass23,weight);
      h_mass13_mc->Fill(mass13,weight);
      if(abs(mode)>=11 && abs(mode)<=13){//CC single pi
        h_mass12_cc_single_pi_mc->Fill(mass12,weight);
        h_mass23_cc_single_pi_mc->Fill(mass23,weight);
        h_mass13_cc_single_pi_mc->Fill(mass13,weight);
      }
      else if(abs(mode)==21){//CC multi pi
        h_mass12_cc_multi_pi_mc->Fill(mass12,weight);
        h_mass23_cc_multi_pi_mc->Fill(mass23,weight);
        h_mass13_cc_multi_pi_mc->Fill(mass13,weight);
      }
      else if(abs(mode)==22){//CC eta production
        h_mass12_cc_eta_mc->Fill(mass12,weight);
        h_mass23_cc_eta_mc->Fill(mass23,weight);
        h_mass13_cc_eta_mc->Fill(mass13,weight);
      }
      else if(abs(mode)==26){//CC DIS production
        h_mass12_cc_dis_mc->Fill(mass12,weight);
        h_mass23_cc_dis_mc->Fill(mass23,weight);
        h_mass13_cc_dis_mc->Fill(mass13,weight);
      }
      else if(abs(mode)==41){//NC multi pi
        h_mass12_nc_multi_pi_mc->Fill(mass12,weight);
        h_mass23_nc_multi_pi_mc->Fill(mass23,weight);
        h_mass13_nc_multi_pi_mc->Fill(mass13,weight);
      }
      else if(abs(mode)==46){//NC DIS production
        h_mass12_nc_dis_mc->Fill(mass12,weight);
        h_mass23_nc_dis_mc->Fill(mass23,weight);
        h_mass13_nc_dis_mc->Fill(mass13,weight);
      }
      else {//Other mode 
        h_mass12_other_mc->Fill(mass12,weight);
        h_mass23_other_mc->Fill(mass23,weight);
        h_mass13_other_mc->Fill(mass13,weight);
      }
      evt_mc_3bin += weight;
    }
    if(total_mass_mc>800 && total_mass_mc<1050) evt_mc_box += weight;

  }
  output->Write();
  cout << "sk" << period << endl;
  float diff_3bin = evt_dt_3bin - evt_mc_3bin;
  float sigma_3bin = diff_3bin/sqrt(evt_mc_3bin);
  cout << "events in 750<M<900MeV data/mc/diff/sigma=" 
    << evt_dt_3bin << "/" << evt_mc_3bin << "/" <<  diff_3bin << "/" << sigma_3bin << endl;
  float diff_box = evt_dt_box - evt_mc_box;
  float sigma_box = diff_box/sqrt(evt_mc_box);
  cout << "events in 800<M<1050MeV data/mc/diff/sigma=" 
    << evt_dt_box << "/" << evt_mc_box << "/" <<  diff_box << "/" << sigma_box << endl;
}



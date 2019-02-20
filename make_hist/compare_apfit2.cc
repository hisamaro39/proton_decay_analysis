void compare_apfit2(){

  TFile *input_miura = TFile::Open("../output/p_mumumu_miura.sk4.mode_p_mumumu_miura_tree.root");
  TTree *tree_miura = (TTree*) input_miura->Get("osc_tuple");
  TFile *input_take = TFile::Open("../output/p_mumumu_take.sk4.mode_p_mumumu_take_tree.root");
  TTree *tree_take = (TTree*) input_take->Get("osc_tuple");
  int nring_miura,nring_take,nmulike_miura,nmulike_take;
  float total_mass_miura,total_mass_take,total_mom_miura,total_mom_take,dlfct_miura,dlfct_take;
  float vertex_x_miura,vertex_x_take;
  float vertex_y_miura,vertex_y_take;
  float vertex_z_miura,vertex_z_take;
  float mmom_min_miura,mmom_min_take,mmom_mid_miura,mmom_mid_take,mmom_max_miura,mmom_max_take;
  float prob_angle_miura[5],prob_angle_take[5],mmom_miura[5],mmom_take[5],prob_pattern_miura[5],prob_pattern_take[5];
  float probms_e_miura[5],probms_e_take[5],probms_mu_miura[5],probms_mu_take[5];
  float prmslg_e_miura[5],prmslg_e_take[5],prmslg_mu_miura[5],prmslg_mu_take[5];
  float ang_miura[5],ang_take[5],ange_miura[5],ange_take[5],angm_miura[5],angm_take[5];
  float dir_x_miura[5],dir_x_take[5];
  float dir_y_miura[5],dir_y_take[5];
  float dir_z_miura[5],dir_z_take[5];
  tree_miura->SetBranchAddress("nring",&nring_miura);
  tree_take->SetBranchAddress("nring",&nring_take);
  tree_miura->SetBranchAddress("nmulike",&nmulike_miura);
  tree_take->SetBranchAddress("nmulike",&nmulike_take);
  tree_miura->SetBranchAddress("total_mass",&total_mass_miura);
  tree_take->SetBranchAddress("total_mass",&total_mass_take);
  tree_miura->SetBranchAddress("total_mom",&total_mom_miura);
  tree_take->SetBranchAddress("total_mom",&total_mom_take);
  tree_miura->SetBranchAddress("dlfct",&dlfct_miura);
  tree_take->SetBranchAddress("dlfct",&dlfct_take);
  tree_miura->SetBranchAddress("prob_angle",&prob_angle_miura);
  tree_take->SetBranchAddress("prob_angle",&prob_angle_take);
  tree_miura->SetBranchAddress("probms_e",&probms_e_miura);
  tree_take->SetBranchAddress("probms_e",&probms_e_take);
  tree_miura->SetBranchAddress("probms_mu",&probms_mu_miura);
  tree_take->SetBranchAddress("probms_mu",&probms_mu_take);
  tree_miura->SetBranchAddress("prob_pattern",&prob_pattern_miura);
  tree_take->SetBranchAddress("prob_pattern",&prob_pattern_take);
  tree_miura->SetBranchAddress("prmslg_e",&prmslg_e_miura);
  tree_take->SetBranchAddress("prmslg_e",&prmslg_e_take);
  tree_miura->SetBranchAddress("prmslg_mu",&prmslg_mu_miura);
  tree_take->SetBranchAddress("prmslg_mu",&prmslg_mu_take);
  tree_miura->SetBranchAddress("dlfct",&dlfct_miura);
  tree_take->SetBranchAddress("dlfct",&dlfct_take);
  tree_miura->SetBranchAddress("mmom",&mmom_miura);
  tree_take->SetBranchAddress("mmom",&mmom_take);
  tree_miura->SetBranchAddress("mmom_min",&mmom_min_miura);
  tree_take->SetBranchAddress("mmom_min",&mmom_min_take);
  tree_miura->SetBranchAddress("mmom_mid",&mmom_mid_miura);
  tree_take->SetBranchAddress("mmom_mid",&mmom_mid_take);
  tree_miura->SetBranchAddress("mmom_max",&mmom_max_miura);
  tree_take->SetBranchAddress("mmom_max",&mmom_max_take);
  tree_miura->SetBranchAddress("vertex_x",&vertex_x_miura);
  tree_take->SetBranchAddress("vertex_x",&vertex_x_take);
  tree_miura->SetBranchAddress("vertex_y",&vertex_y_miura);
  tree_take->SetBranchAddress("vertex_y",&vertex_y_take);
  tree_miura->SetBranchAddress("vertex_z",&vertex_z_miura);
  tree_take->SetBranchAddress("vertex_z",&vertex_z_take);
  tree_miura->SetBranchAddress("dir_x",&dir_x_miura);
  tree_take->SetBranchAddress("dir_x",&dir_x_take);
  tree_miura->SetBranchAddress("dir_y",&dir_y_miura);
  tree_take->SetBranchAddress("dir_y",&dir_y_take);
  tree_miura->SetBranchAddress("dir_z",&dir_z_miura);
  tree_take->SetBranchAddress("dir_z",&dir_z_take);
  tree_miura->SetBranchAddress("ang",&ang_miura);
  tree_take->SetBranchAddress("ang",&ang_take);
  tree_miura->SetBranchAddress("angm",&angm_miura);
  tree_take->SetBranchAddress("angm",&angm_take);

  TFile *output = new TFile("output/output_compare_apfit_miura_take_p_mumumu.root","recreate");
  TH1 *h_diff_dlfct = new TH1F("h_diff_dlfct",";;",100,-1,1);
  TH1 *h_diff_dlfct_same_vertex = new TH1F("h_diff_dlfct_same_vertex",";;",100,-1,1);
  TH1 *h_diff_dlfct_different_vertex = new TH1F("h_diff_dlfct_different_vertex",";;",100,-1,1);
  TH1 *h_diff_vertex_x = new TH1F("h_diff_vertex_x",";;",100,-100,100);
  TH1 *h_diff_vertex_y = new TH1F("h_diff_vertex_y",";;",100,-100,100);
  TH1 *h_diff_vertex_z = new TH1F("h_diff_vertex_z",";;",100,-100,100);
  TH1 *h_diff_total_mass = new TH1F("h_diff_total_mass",";;",100,-100,100);
  TH1 *h_diff_total_mom = new TH1F("h_diff_total_mom",";;",100,-100,100);
  TH1 *h_diff_mmom = new TH1F("h_diff_mmom",";;",100,-100,100);
  TH1 *h_diff_mmom_same_vertex = new TH1F("h_diff_mmom_same_vertex",";;",100,-100,100);
  TH1 *h_diff_mmom_different_vertex = new TH1F("h_diff_mmom_different_vertex",";;",100,-100,100);
  TH1 *h_diff_dir_x = new TH1F("h_diff_dir_x",";;",100,-1,1);
  TH1 *h_diff_dir_y = new TH1F("h_diff_dir_y",";;",100,-1,1);
  TH1 *h_diff_dir_z = new TH1F("h_diff_dir_z",";;",100,-1,1);
  TH1 *h_diff_dir_x_same_vertex = new TH1F("h_diff_dir_x_same_vertex",";;",100,-1,1);
  TH1 *h_diff_dir_y_same_vertex = new TH1F("h_diff_dir_y_same_vertex",";;",100,-1,1);
  TH1 *h_diff_dir_z_same_vertex = new TH1F("h_diff_dir_z_same_vertex",";;",100,-1,1);
  TH1 *h_diff_prob_angle = new TH1F("h_diff_prob_angle",";;",100,-1,1);
  TH1 *h_diff_prob_angle_same_vertex = new TH1F("h_diff_prob_angle_same_vertex",";;",100,-1,1);
  TH1 *h_diff_prob_angle_different_vertex = new TH1F("h_diff_prob_angle_different_vertex",";;",100,-1,1);
  TH1 *h_diff_probms_e = new TH1F("h_diff_probms_e",";;",100,-1,1);
  TH1 *h_diff_probms_mu = new TH1F("h_diff_probms_mu",";;",100,-1,1);
  TH1 *h_diff_prob_pattern = new TH1F("h_diff_prob_pattern",";;",100,-100,100);
  TH1 *h_diff_prmslg_e = new TH1F("h_diff_prmslg_e",";;",100,-100,100);
  TH1 *h_diff_prmslg_mu = new TH1F("h_diff_prmslg_mu",";;",100,-100,100);
  TH1 *h_diff_mmom_min = new TH1F("h_diff_mmom_min",";;",100,-100,100);
  TH1 *h_diff_mmom_mid = new TH1F("h_diff_mmom_mid",";;",100,-100,100);
  TH1 *h_diff_mmom_max = new TH1F("h_diff_mmom_max",";;",100,-100,100);
  TH1 *h_diff_ang = new TH1F("h_diff_ang",";;",100,-20,20);
  TH1 *h_diff_ang_same_vertex = new TH1F("h_diff_ang_same_vertex",";;",100,-20,20);
  TH1 *h_diff_ang_different_vertex = new TH1F("h_diff_ang_different_vertex",";;",100,-20,20);
  TH1 *h_diff_angm = new TH1F("h_diff_angm",";;",100,-20,20);
  TH1 *h_diff_angm_same_vertex = new TH1F("h_diff_angm_same_vertex",";;",100,-20,20);
  TH1 *h_diff_angm_different_vertex = new TH1F("h_diff_angm_different_vertex",";;",100,-20,20);

  for (int e=0;e<tree_miura->GetEntries();e++){
    cout << "event" << e << endl;
    tree_miura->GetEntry(e);
    tree_take->GetEntry(e);
    cout << "nring miura/take=" << nring_miura << "/" << nring_take << endl;
    cout << "nmulike miura/take=" << nmulike_miura << "/" << nmulike_take << endl;
    cout << "total_mass miura/take=" << total_mass_miura << "/" << total_mass_take << endl;
    cout << "total_mom miura/take=" << total_mom_miura << "/" << total_mom_take << endl;
    cout << "dlfct miura/take=" << dlfct_miura << "/" << dlfct_take << endl;
    cout << "vertex_x miura/take=" << vertex_x_miura << "/" << vertex_x_take << endl;
    cout << "vertex_y miura/take=" << vertex_y_miura << "/" << vertex_y_take << endl;
    cout << "vertex_z miura/take=" << vertex_z_miura << "/" << vertex_z_take << endl;
    float diff_dlfct = dlfct_take - dlfct_miura;
    float diff_vertex_x = vertex_x_take - vertex_x_miura;
    float diff_vertex_y = vertex_y_take - vertex_y_miura;
    float diff_vertex_z = vertex_z_take - vertex_z_miura;
    h_diff_dlfct->Fill(diff_dlfct);
    cout << "diff_vertex x/y/z=" << diff_vertex_x << "/" << diff_vertex_y << "/" << diff_vertex_z << endl;
    if(fabs(diff_vertex_x) < 1e-5 && fabs(diff_vertex_y) < 1e-5 && fabs(diff_vertex_z) < 1e-5) {
      h_diff_dlfct_same_vertex->Fill(diff_dlfct);
    }
    else {
      h_diff_dlfct_different_vertex->Fill(diff_dlfct);
    }
    h_diff_vertex_x->Fill(diff_vertex_x);
    h_diff_vertex_y->Fill(diff_vertex_y);
    h_diff_vertex_z->Fill(diff_vertex_z);
    if(nring_miura==nring_take && nmulike_miura==nmulike_take){
      if(nring_miura==3 && nmulike_miura==3){
        h_diff_total_mass->Fill(total_mass_take-total_mass_miura);
        h_diff_total_mom->Fill(total_mom_take-total_mom_miura);
        h_diff_mmom_min->Fill(mmom_min_take-mmom_min_miura);
        h_diff_mmom_mid->Fill(mmom_mid_take-mmom_mid_miura);
        h_diff_mmom_max->Fill(mmom_max_take-mmom_max_miura);
      }
      for(int r=0;r<nring_miura;r++){
        cout << "ring" << r << endl;
        cout << "prob_angle miura/take=" << prob_angle_miura[r] << "/" << prob_angle_take[r] << endl;
        cout << "probms_e miura/take=" << probms_e_miura[r] << "/" << probms_e_take[r] << endl;
        cout << "probms_mu miura/take=" << probms_mu_miura[r] << "/" << probms_mu_take[r] << endl;
        cout << "prob_pattern miura/take=" << prob_pattern_miura[r] << "/" << prob_pattern_take[r] << endl;
        cout << "prmslg_e miura/take=" << prmslg_e_miura[r] << "/" << prmslg_e_take[r] << endl;
        cout << "prmslg_mu miura/take=" << prmslg_mu_miura[r] << "/" << prmslg_mu_take[r] << endl;
        cout << "mmom miura/take=" << mmom_miura[r] << "/" << mmom_take[r] << endl;
        cout << "dir_x miura/take=" << dir_x_miura[r] << "/" << dir_x_take[r] << endl;
        cout << "dir_y miura/take=" << dir_y_miura[r] << "/" << dir_y_take[r] << endl;
        cout << "dir_z miura/take=" << dir_z_miura[r] << "/" << dir_z_take[r] << endl;
        cout << "ang miura/take=" << ang_miura[r] << "/" << ang_take[r] << endl;
        cout << "angm miura/take=" << angm_miura[r] << "/" << angm_take[r] << endl;
        float diff_mmom = mmom_take[r]-mmom_miura[r];
        float diff_prob_angle = prob_angle_take[r]-prob_angle_miura[r];
        float diff_probms_e = sqrt(fabs(probms_e_take[r]))-sqrt(fabs(probms_e_miura[r]));
        float diff_probms_mu = sqrt(fabs(probms_mu_take[r]))-sqrt(fabs(probms_mu_miura[r]));
        float diff_prob_pattern = prob_pattern_take[r]-prob_pattern_miura[r];
        float diff_prmslg_e = sqrt(fabs(prmslg_e_take[r]))-sqrt(fabs(prmslg_e_miura[r]));
        float diff_prmslg_mu = sqrt(fabs(prmslg_mu_take[r]))-sqrt(fabs(prmslg_mu_miura[r]));
        float diff_dir_x = dir_x_take[r]-dir_x_miura[r];
        float diff_dir_y = dir_y_take[r]-dir_y_miura[r];
        float diff_dir_z = dir_z_take[r]-dir_z_miura[r];
        float diff_ang = ang_take[r]-ang_miura[r];
        float diff_angm = angm_take[r]-angm_miura[r];
        h_diff_mmom->Fill(diff_mmom);
        h_diff_prob_angle->Fill(diff_prob_angle);
        h_diff_probms_e->Fill(diff_probms_e);
        h_diff_probms_mu->Fill(diff_probms_mu);
        h_diff_prob_pattern->Fill(diff_prob_pattern);
        h_diff_prmslg_e->Fill(diff_prmslg_e);
        h_diff_prmslg_mu->Fill(diff_prmslg_mu);
        h_diff_dir_x->Fill(diff_dir_x);
        h_diff_dir_y->Fill(diff_dir_y);
        h_diff_dir_z->Fill(diff_dir_z);
        h_diff_ang->Fill(diff_ang);
        h_diff_angm->Fill(diff_angm);
        if(fabs(diff_vertex_x) < 1e-5 && fabs(diff_vertex_y) < 1e-5 && fabs(diff_vertex_z) < 1e-5){
          h_diff_dir_x_same_vertex->Fill(diff_dir_x);
          h_diff_dir_y_same_vertex->Fill(diff_dir_y);
          h_diff_dir_z_same_vertex->Fill(diff_dir_z);
          h_diff_prob_angle_same_vertex->Fill(prob_angle_take[r]-prob_angle_miura[r]);
          h_diff_ang_same_vertex->Fill(diff_ang);
          h_diff_angm_same_vertex->Fill(diff_angm);
          h_diff_mmom_same_vertex->Fill(diff_mmom);
        }else {
          h_diff_prob_angle_different_vertex->Fill(prob_angle_take[r]-prob_angle_miura[r]);
          h_diff_ang_different_vertex->Fill(diff_ang);
          h_diff_angm_different_vertex->Fill(diff_angm);
          h_diff_mmom_different_vertex->Fill(diff_mmom);
        }
      }
    }
  }
  output->Write();

}



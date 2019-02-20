void compare_apfit(){

  TFile *input_def = TFile::Open("../output/p_mumumu_def.sk4.mode_p_mumumu_def_tree.root");
  TTree *tree_def = (TTree*) input_def->Get("osc_tuple");
  TFile *input_miura = TFile::Open("../output/p_mumumu_miura.sk4.mode_p_mumumu_miura_tree.root");
  TTree *tree_miura = (TTree*) input_miura->Get("osc_tuple");
  int nring_def,nring_miura,nmulike_def,nmulike_miura;
  float total_mass_def,total_mass_miura,total_mom_def,total_mom_miura,dlfct_def,dlfct_miura;
  float vertex_x_def,vertex_x_miura;
  float vertex_y_def,vertex_y_miura;
  float vertex_z_def,vertex_z_miura;
  float mmom_min_def,mmom_min_miura,mmom_mid_def,mmom_mid_miura,mmom_max_def,mmom_max_miura;
  float prob_angle_def[5],prob_angle_miura[5],mmom_def[5],mmom_miura[5],prob_pattern_def[5],prob_pattern_miura[5];
  float probms_e_def[5],probms_e_miura[5],probms_mu_def[5],probms_mu_miura[5];
  float prmslg_e_def[5],prmslg_e_miura[5],prmslg_mu_def[5],prmslg_mu_miura[5];
  float ang_def[5],ang_miura[5],ange_def[5],ange_miura[5],angm_def[5],angm_miura[5];
  float dir_x_def[5],dir_x_miura[5];
  float dir_y_def[5],dir_y_miura[5];
  float dir_z_def[5],dir_z_miura[5];
  tree_def->SetBranchAddress("nring",&nring_def);
  tree_miura->SetBranchAddress("nring",&nring_miura);
  tree_def->SetBranchAddress("nmulike",&nmulike_def);
  tree_miura->SetBranchAddress("nmulike",&nmulike_miura);
  tree_def->SetBranchAddress("total_mass",&total_mass_def);
  tree_miura->SetBranchAddress("total_mass",&total_mass_miura);
  tree_def->SetBranchAddress("total_mom",&total_mom_def);
  tree_miura->SetBranchAddress("total_mom",&total_mom_miura);
  tree_def->SetBranchAddress("dlfct",&dlfct_def);
  tree_miura->SetBranchAddress("dlfct",&dlfct_miura);
  tree_def->SetBranchAddress("prob_angle",&prob_angle_def);
  tree_miura->SetBranchAddress("prob_angle",&prob_angle_miura);
  tree_def->SetBranchAddress("probms_e",&probms_e_def);
  tree_miura->SetBranchAddress("probms_e",&probms_e_miura);
  tree_def->SetBranchAddress("probms_mu",&probms_mu_def);
  tree_miura->SetBranchAddress("probms_mu",&probms_mu_miura);
  tree_def->SetBranchAddress("prob_pattern",&prob_pattern_def);
  tree_miura->SetBranchAddress("prob_pattern",&prob_pattern_miura);
  tree_def->SetBranchAddress("prmslg_e",&prmslg_e_def);
  tree_miura->SetBranchAddress("prmslg_e",&prmslg_e_miura);
  tree_def->SetBranchAddress("prmslg_mu",&prmslg_mu_def);
  tree_miura->SetBranchAddress("prmslg_mu",&prmslg_mu_miura);
  tree_def->SetBranchAddress("dlfct",&dlfct_def);
  tree_miura->SetBranchAddress("dlfct",&dlfct_miura);
  tree_def->SetBranchAddress("mmom",&mmom_def);
  tree_miura->SetBranchAddress("mmom",&mmom_miura);
  tree_def->SetBranchAddress("mmom_min",&mmom_min_def);
  tree_miura->SetBranchAddress("mmom_min",&mmom_min_miura);
  tree_def->SetBranchAddress("mmom_mid",&mmom_mid_def);
  tree_miura->SetBranchAddress("mmom_mid",&mmom_mid_miura);
  tree_def->SetBranchAddress("mmom_max",&mmom_max_def);
  tree_miura->SetBranchAddress("mmom_max",&mmom_max_miura);
  tree_def->SetBranchAddress("vertex_x",&vertex_x_def);
  tree_miura->SetBranchAddress("vertex_x",&vertex_x_miura);
  tree_def->SetBranchAddress("vertex_y",&vertex_y_def);
  tree_miura->SetBranchAddress("vertex_y",&vertex_y_miura);
  tree_def->SetBranchAddress("vertex_z",&vertex_z_def);
  tree_miura->SetBranchAddress("vertex_z",&vertex_z_miura);
  tree_def->SetBranchAddress("dir_x",&dir_x_def);
  tree_miura->SetBranchAddress("dir_x",&dir_x_miura);
  tree_def->SetBranchAddress("dir_y",&dir_y_def);
  tree_miura->SetBranchAddress("dir_y",&dir_y_miura);
  tree_def->SetBranchAddress("dir_z",&dir_z_def);
  tree_miura->SetBranchAddress("dir_z",&dir_z_miura);
  tree_def->SetBranchAddress("ang",&ang_def);
  tree_miura->SetBranchAddress("ang",&ang_miura);
  tree_def->SetBranchAddress("angm",&angm_def);
  tree_miura->SetBranchAddress("angm",&angm_miura);

  TFile *output = new TFile("output/output_compare_apfit_def_miura_p_mumumu.root","recreate");
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

  for (int e=0;e<tree_def->GetEntries();e++){
    cout << "event" << e << endl;
    tree_def->GetEntry(e);
    tree_miura->GetEntry(e);
    cout << "nring def/miura=" << nring_def << "/" << nring_miura << endl;
    cout << "nmulike def/miura=" << nmulike_def << "/" << nmulike_miura << endl;
    cout << "total_mass def/miura=" << total_mass_def << "/" << total_mass_miura << endl;
    cout << "total_mom def/miura=" << total_mom_def << "/" << total_mom_miura << endl;
    cout << "dlfct def/miura=" << dlfct_def << "/" << dlfct_miura << endl;
    cout << "vertex_x def/miura=" << vertex_x_def << "/" << vertex_x_miura << endl;
    cout << "vertex_y def/miura=" << vertex_y_def << "/" << vertex_y_miura << endl;
    cout << "vertex_z def/miura=" << vertex_z_def << "/" << vertex_z_miura << endl;
    float diff_dlfct = dlfct_miura - dlfct_def;
    float diff_vertex_x = vertex_x_miura - vertex_x_def;
    float diff_vertex_y = vertex_y_miura - vertex_y_def;
    float diff_vertex_z = vertex_z_miura - vertex_z_def;
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
    if(nring_def==nring_miura && nmulike_def==nmulike_miura){
      if(nring_def==3 && nmulike_def==3){
        h_diff_total_mass->Fill(total_mass_miura-total_mass_def);
        h_diff_total_mom->Fill(total_mom_miura-total_mom_def);
        h_diff_mmom_min->Fill(mmom_min_miura-mmom_min_def);
        h_diff_mmom_mid->Fill(mmom_mid_miura-mmom_mid_def);
        h_diff_mmom_max->Fill(mmom_max_miura-mmom_max_def);
      }
      for(int r=0;r<nring_def;r++){
        cout << "ring" << r << endl;
        cout << "prob_angle def/miura=" << prob_angle_def[r] << "/" << prob_angle_miura[r] << endl;
        cout << "probms_e def/miura=" << probms_e_def[r] << "/" << probms_e_miura[r] << endl;
        cout << "probms_mu def/miura=" << probms_mu_def[r] << "/" << probms_mu_miura[r] << endl;
        cout << "prob_pattern def/miura=" << prob_pattern_def[r] << "/" << prob_pattern_miura[r] << endl;
        cout << "prmslg_e def/miura=" << prmslg_e_def[r] << "/" << prmslg_e_miura[r] << endl;
        cout << "prmslg_mu def/miura=" << prmslg_mu_def[r] << "/" << prmslg_mu_miura[r] << endl;
        cout << "mmom def/miura=" << mmom_def[r] << "/" << mmom_miura[r] << endl;
        cout << "dir_x def/miura=" << dir_x_def[r] << "/" << dir_x_miura[r] << endl;
        cout << "dir_y def/miura=" << dir_y_def[r] << "/" << dir_y_miura[r] << endl;
        cout << "dir_z def/miura=" << dir_z_def[r] << "/" << dir_z_miura[r] << endl;
        cout << "ang def/miura=" << ang_def[r] << "/" << ang_miura[r] << endl;
        cout << "angm def/miura=" << angm_def[r] << "/" << angm_miura[r] << endl;
        float diff_mmom = mmom_miura[r]-mmom_def[r];
        float diff_prob_angle = prob_angle_miura[r]-prob_angle_def[r];
        float diff_probms_e = sqrt(fabs(probms_e_miura[r]))-sqrt(fabs(probms_e_def[r]));
        float diff_probms_mu = sqrt(fabs(probms_mu_miura[r]))-sqrt(fabs(probms_mu_def[r]));
        float diff_prob_pattern = prob_pattern_miura[r]-prob_pattern_def[r];
        float diff_prmslg_e = sqrt(fabs(prmslg_e_miura[r]))-sqrt(fabs(prmslg_e_def[r]));
        float diff_prmslg_mu = sqrt(fabs(prmslg_mu_miura[r]))-sqrt(fabs(prmslg_mu_def[r]));
        float diff_dir_x = dir_x_miura[r]-dir_x_def[r];
        float diff_dir_y = dir_y_miura[r]-dir_y_def[r];
        float diff_dir_z = dir_z_miura[r]-dir_z_def[r];
        float diff_ang = ang_miura[r]-ang_def[r];
        float diff_angm = angm_miura[r]-angm_def[r];
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
          h_diff_prob_angle_same_vertex->Fill(prob_angle_miura[r]-prob_angle_def[r]);
          h_diff_ang_same_vertex->Fill(diff_ang);
          h_diff_angm_same_vertex->Fill(diff_angm);
          h_diff_mmom_same_vertex->Fill(diff_mmom);
        }else {
          h_diff_prob_angle_different_vertex->Fill(prob_angle_miura[r]-prob_angle_def[r]);
          h_diff_ang_different_vertex->Fill(diff_ang);
          h_diff_angm_different_vertex->Fill(diff_angm);
          h_diff_mmom_different_vertex->Fill(diff_mmom);
        }
      }
    }
  }
  output->Write();

}



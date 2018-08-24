void make_efficiency_plot(){

  gROOT->SetStyle("Plain");
  gStyle->SetOptStat(0);
  
  TFile *input[4];
  input[0] = TFile::Open("../output/p_epi.sk4.18a.root");
  input[1] = TFile::Open("../output/p_mupi.sk4.18a.root");
  input[2] = TFile::Open("../output/p_eee.sk4.18a.root");
  input[3] = TFile::Open("../output/p_mumumu.sk4.18a.root");
  TFile *output = new TFile("output/efficiency.root","recreate");
  vector<string> bunbo,bunsi;
  vector<int> file_type;
  bunbo.clear();bunsi.clear();file_type.clear();
  bunbo.push_back("mom_gamma_vector_only_epi0");bunsi.push_back("mom_gamma_match_ring_vector_only_epi0");
  file_type.push_back(0);
  bunbo.push_back("mom_gamma_match_ring_vector_only_epi0");bunsi.push_back("mom_gamma_vector_prob_ring_angle_elike_gamma_match_ring_only_epi0");
  file_type.push_back(0);
  bunbo.push_back("mom_gamma_match_ring_vector_only_epi0");bunsi.push_back("mom_gamma_vector_prob_hit_elike_gamma_match_ring_only_epi0");
  file_type.push_back(0);
  bunbo.push_back("mom_gamma_match_ring_vector_only_epi0");bunsi.push_back("mom_gamma_vector_pid_elike_gamma_match_ring_only_epi0");
  file_type.push_back(0);
  bunbo.push_back("mom_e_match_ring_vector_only_epi0");bunsi.push_back("mom_e_vector_prob_ring_angle_elike_e_match_ring_only_epi0");
  file_type.push_back(0);
  bunbo.push_back("mom_e_match_ring_vector_only_epi0");bunsi.push_back("mom_e_vector_prob_hit_elike_e_match_ring_only_epi0");
  file_type.push_back(0);
  bunbo.push_back("mom_mu_match_ring_vector_only_mupi0");bunsi.push_back("mom_mu_vector_prob_ring_angle_mulike_mu_match_ring_only_mupi0");
  file_type.push_back(1);
  bunbo.push_back("mom_mu_match_ring_vector_only_mupi0");bunsi.push_back("mom_mu_vector_prob_hit_mulike_mu_match_ring_only_mupi0");
  file_type.push_back(1);
  bunbo.push_back("mom_e_match_ring_vector_only_eee");bunsi.push_back("mom_e_vector_prob_ring_angle_elike_e_match_ring_only_eee");
  file_type.push_back(2);
  bunbo.push_back("mom_e_match_ring_vector_only_eee");bunsi.push_back("mom_e_vector_prob_hit_elike_e_match_ring_only_eee");
  file_type.push_back(2);
  bunbo.push_back("mom_e_match_ring_vector_only_eee");bunsi.push_back("mom_e_vector_pid_elike_e_match_ring_only_eee");
  file_type.push_back(2);
  bunbo.push_back("mom_mu_match_ring_vector_only_mumumu");bunsi.push_back("mom_mu_vector_prob_ring_angle_mulike_mu_match_ring_only_mumumu");
  file_type.push_back(3);
  bunbo.push_back("mom_mu_match_ring_vector_only_mumumu");bunsi.push_back("mom_mu_vector_prob_hit_mulike_mu_match_ring_only_mumumu");
  file_type.push_back(3);
  bunbo.push_back("mom_mu_match_ring_vector_only_mumumu");bunsi.push_back("mom_mu_vector_pid_mulike_mu_match_ring_only_mumumu");
  file_type.push_back(3);
  bunbo.push_back("mom_e_vector_only_eee");bunsi.push_back("mom_e_match_ring_vector_only_eee");
  file_type.push_back(2);
  bunbo.push_back("mom_mu_vector_only_mumumu");bunsi.push_back("mom_mu_match_ring_vector_only_mumumu");
  file_type.push_back(3);

  for(int n=0;n<bunbo.size();n++){
    cout << "bunbo/bunsi=" << bunbo[n] << "/" << bunsi[n] << endl;
    TH1 *h_bunbo = (TH1*) input[file_type[n]]->Get(bunbo[n].c_str());
    TH1 *h_bunsi = (TH1*) input[file_type[n]]->Get(bunsi[n].c_str());
    TGraphAsymmErrors *eff = new TGraphAsymmErrors();
    eff->Divide(h_bunsi,h_bunbo);
    eff->SetName(Form("efficiency_%s",bunsi[n].c_str()));
    //h_bunbo->Draw();
    //h_bunsi->SetLineColor(2);
    //h_bunsi->Draw("same");
    //eff->Draw("ap");
    cout << eff->GetName() << endl;
    output->Add(eff);
  }
  output->Write();

}

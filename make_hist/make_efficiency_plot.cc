void make_efficiency_plot(){

  string type[] = {"p_epi","p_mupi","p_eee","p_mumumu","p_emumu","p_muee","fcmc","fcdt"};
  gROOT->SetStyle("Plain");
  gStyle->SetOptStat(0);
  
  TFile *output = new TFile("output/efficiency.root","recreate");
  vector<string> bunbo,bunsi;
  vector<int> input_type,mode_type;
  input_type.clear();mode_type.clear();
  bunbo.clear();bunsi.clear();
  bunbo.push_back("true_mom_lepton");bunsi.push_back("true_mom_lepton_match_ring");
  input_type.push_back(2);mode_type.push_back(2);
  bunbo.push_back("true_mom_lepton");bunsi.push_back("true_mom_lepton_match_ring");
  input_type.push_back(3);mode_type.push_back(3);
  bunbo.push_back("true_mom_lepton_nring3");bunsi.push_back("true_mom_lepton_match_ring_nring3");
  input_type.push_back(2);mode_type.push_back(2);
  bunbo.push_back("true_mom_lepton_nring3");bunsi.push_back("true_mom_lepton_match_ring_nring3");
  input_type.push_back(3);mode_type.push_back(3);
  bunbo.push_back("true_mom_muon_nring3");bunsi.push_back("true_mom_muon_match_ring_nring3");
  input_type.push_back(5);mode_type.push_back(5);
  bunbo.push_back("true_mom_lepton_match_ring");bunsi.push_back("true_mom_lepton_match_ring_angle_mulike");
  input_type.push_back(3);mode_type.push_back(3);
  bunbo.push_back("true_mom_lepton_match_ring");bunsi.push_back("true_mom_lepton_match_ring_charge_mulike");
  input_type.push_back(3);mode_type.push_back(3);
  bunbo.push_back("true_mom_lepton_match_ring");bunsi.push_back("true_mom_lepton_match_ring_angle_elike");
  input_type.push_back(2);mode_type.push_back(2);
  bunbo.push_back("true_mom_lepton_match_ring");bunsi.push_back("true_mom_lepton_match_ring_charge_elike");
  input_type.push_back(2);mode_type.push_back(2);
  bunbo.push_back("true_mom_lepton_match_ring_nring3");bunsi.push_back("true_mom_lepton_match_ring_angle_elike_nring3");
  input_type.push_back(2);mode_type.push_back(2);
  bunbo.push_back("true_mom_lepton_match_ring_nring3");bunsi.push_back("true_mom_lepton_match_ring_angle_mulike_nring3");
  input_type.push_back(3);mode_type.push_back(3);
  bunbo.push_back("true_mom_muon_match_ring_nring3");bunsi.push_back("true_mom_muon_match_ring_angle_mulike_nring3");
  input_type.push_back(5);mode_type.push_back(5);
  bunbo.push_back("true_mom_muon_match_ring_nring3");bunsi.push_back("true_mom_muon_match_ring_angle_mulike_nring3");
  input_type.push_back(1);mode_type.push_back(1);

  TFile *input;
  for(int n=0;n<bunbo.size();n++){
    //cout << "bunbo/bunsi=" << bunbo[n] << "/" << bunsi[n] << endl;
    input = TFile::Open(Form("../output/%s.sk4.mode_%s_validation.root",type[input_type[n]].c_str(),type[mode_type[n]].c_str()));
    TGraphAsymmErrors *eff = new TGraphAsymmErrors();
    eff->SetName(Form("efficiency_%s_%s",bunsi[n].c_str(),type[mode_type[n]].c_str()));
    TH1 *h_bunbo = (TH1*) input->Get(bunbo[n].c_str());
    TH1 *h_bunsi = (TH1*) input->Get(bunsi[n].c_str());
    int nbins = h_bunbo->GetNbinsX();
    for(int b=0;b<nbins;b++){
      float bin_center = h_bunbo->GetBinCenter(b+1);
      float bin_width = h_bunbo->GetBinWidth(b+1);
      //cout << "bin center/width=" << bin_center << "/" << bin_width << endl;
      float n_bunsi = h_bunsi->GetBinContent(b+1);
      float n_bunbo = h_bunbo->GetBinContent(b+1);
      float ratio = (n_bunbo)? n_bunsi/n_bunbo : 0;
      float err = (n_bunbo)? sqrt(n_bunsi*n_bunbo*(n_bunsi+n_bunbo))/pow(n_bunbo,2) : 0;
      float err_high = (ratio+err>1)? 1-ratio : err;
      float err_low = (ratio>err)? err : ratio; 
      //cout << "bunsi/bunbo/ratio=" << n_bunsi << "/" << n_bunbo << "/" << ratio << "+-" << err << endl;
      eff->SetPoint(b+1,bin_center,ratio);
      eff->SetPointError(b+1,bin_width*0.5,bin_width*0.5,err_low,err_high);
    }
    //eff->Divide(h_bunsi,h_bunbo);
    //h_bunbo->Draw();
    //h_bunsi->SetLineColor(2);
    //h_bunsi->Draw("same");
    //eff->Draw("ap");
    //cout << eff->GetName() << endl;
    output->Add(eff);
  }
  output->Write();

}

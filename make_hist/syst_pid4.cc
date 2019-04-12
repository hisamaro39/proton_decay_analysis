void syst_pid4(){
  //setup
  bool calc_chi2 = true;
  int ring_number = 3;

  string str_ring;
  if(ring_number==1) str_ring = "1st";
  if(ring_number==2) str_ring = "2nd";
  if(ring_number==3) str_ring = "3rd";

  TFile *input_mc = TFile::Open("../output/fcmc.sk4.mode_subgev_multiring_tree.root");
  TTree *tree_mc = (TTree*) input_mc->Get("osc_tuple");
  TFile *input_data = TFile::Open("../output/fcdt.sk4.mode_subgev_multiring_tree.root");
  TTree *tree_data = (TTree*) input_data->Get("osc_tuple");

  int nring_data,nring_mc;
  float prob_angle_data[5],prob_angle_mc[5],weight;
  tree_data->SetBranchAddress("nring",&nring_data);
  tree_mc->SetBranchAddress("nring",&nring_mc);
  tree_data->SetBranchAddress("prob_angle",&prob_angle_data);
  tree_mc->SetBranchAddress("prob_angle",&prob_angle_mc);
  tree_mc->SetBranchAddress("weight",&weight);

  TH1 *h_prob_angle_data = new TH1F(Form("h_prob_angle_%s_data",str_ring.c_str()),";;",100,-10,10);
  TH1 *h_prob_angle_mc = new TH1F(Form("h_prob_angle_%s_mc",str_ring.c_str()),";;",100,-10,10);
  TH1 *h_prob_angle_mc_best_fit = new TH1F(Form("h_prob_angle_%s_mc_best_fit",str_ring.c_str()),";;",100,-10,10);

  //data
  for (int e=0;e<tree_data->GetEntries();e++){
    tree_data->GetEntry(e);
    if(nring_data==3) h_prob_angle_data->Fill(prob_angle_data[ring_number-1]);
  }

  //mc
  for (int e=0;e<tree_mc->GetEntries();e++){
    tree_mc->GetEntry(e);
    if(nring_mc==3) h_prob_angle_mc->Fill(prob_angle_mc[ring_number-1],weight);
  }
  float norm = h_prob_angle_data->Integral()/h_prob_angle_mc->Integral();
  cout << "mc normalization=" << norm << endl;

  TCanvas *c = new TCanvas("c1","",800,600);
  h_prob_angle_data->SetLineWidth(2);
  h_prob_angle_data->Draw("E0");
  h_prob_angle_mc->Scale(norm);
  h_prob_angle_mc->SetLineWidth(2);
  h_prob_angle_mc->SetLineColor(2);
  h_prob_angle_mc->Draw("hist same");

  if(calc_chi2){//calculate chi2
    int nbin = h_prob_angle_data->GetNbinsX();
    float min_chi2=9999.,min_scale=0.,min_shift=0.;
    for (int a=0;a<100;a++){
      for(int b=0;b<100;b++){
        TH1 *h_prob_angle_mc_scale_shift = new TH1F(Form("h_prob_angle_%s_mc_scale%d_shift%d",str_ring.c_str(),a,b),";;",100,-10,10);
        h_prob_angle_mc_scale_shift->SetLineColor(4);
        float scale = 0.5 + 0.01*a;
        float shift = 0.01*b -0.5;
        //float scale = 1;
        //float shift = 0;
        for (int e=0;e<tree_mc->GetEntries();e++){
          tree_mc->GetEntry(e);
          if(nring_mc==3) h_prob_angle_mc_scale_shift->Fill(prob_angle_mc[ring_number-1]*scale+shift,weight*norm);
        }
        cout << "scale/shift=" << scale << "/" << shift << endl;
        float chi2=0.;
        for(int bb=0;bb<nbin;bb++){
          float event_data = h_prob_angle_data->GetBinContent(bb+1);
          float event_mc = h_prob_angle_mc_scale_shift->GetBinContent(bb+1);
          float sigma = sqrt(event_data);
          //cout << "event data/mc=" << event_data << "/" << event_mc << endl;
          if(event_data<1E-5) continue;
          float this_diff = pow((event_data-event_mc)/sigma,2);
          //cout << "this_diff=" << this_diff << endl;
          chi2 += this_diff;
        }
        chi2 = chi2/nbin;
        cout << "chi2=" << chi2 << endl;
        if(chi2<min_chi2){
          min_chi2 = chi2;
          min_shift = shift;
          min_scale = scale;
          //h_prob_angle_mc_best_fit = h_prob_angle_mc_scale_shift;
        }
      }
    }
    cout << "min chi2/scale/shift=" << min_chi2 << "/" << min_scale << "/" << min_shift << endl;
  }

  float final_scale,final_shift;
  if(ring_number==1){
    final_scale = 1.49;
    final_shift = -0.49;
  }
  if(ring_number==2){
    final_scale = 0.89;
    final_shift = -0.08;
  }
  if(ring_number==3){
    final_scale = 0.96;
    final_shift = 0.04;
  }
  for (int e=0;e<tree_mc->GetEntries();e++){
    tree_mc->GetEntry(e);
    if(nring_mc==3) h_prob_angle_mc_best_fit->Fill(prob_angle_mc[ring_number-1]*final_scale+final_shift,weight*norm);
  }
  h_prob_angle_mc_best_fit->SetLineColor(4);
  h_prob_angle_mc_best_fit->SetLineWidth(2);
  h_prob_angle_mc_best_fit->Draw("same hist");
  if(!calc_chi2) c->SaveAs(Form("hist/compare_data_mc_prob_angle_before_after_fit_subgev_multiring_%s.pdf",str_ring.c_str()));

}



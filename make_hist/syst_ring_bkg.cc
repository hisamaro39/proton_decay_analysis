void syst_ring_bkg(){
  //setup
  bool calc_chi2 = false;
  int rc_number = 3;
  int period = 4;

  gROOT->SetStyle("Plain");
  gStyle->SetOptStat(0);

  TFile *input_mc = TFile::Open(Form("../output/fcmc_rc.sk%d.mode_subgev_multiring_tree.root",period));
  TTree *tree_mc = (TTree*) input_mc->Get("osc_tuple");
  TFile *input_data = TFile::Open(Form("../output/fcdt_rc.sk%d.mode_subgev_multiring_tree.root",period));
  TTree *tree_data = (TTree*) input_data->Get("osc_tuple");

  int nring_data,nring_mc;
  float rc1_data,rc2_data,rc3_data,rc4_data,rc1_mc,rc2_mc,rc3_mc,rc4_mc,weight;
  tree_data->SetBranchAddress("nring",&nring_data);
  tree_mc->SetBranchAddress("nring",&nring_mc);
  tree_data->SetBranchAddress("rc1",&rc1_data);
  tree_mc->SetBranchAddress("rc1",&rc1_mc);
  tree_data->SetBranchAddress("rc2",&rc2_data);
  tree_mc->SetBranchAddress("rc2",&rc2_mc);
  tree_data->SetBranchAddress("rc3",&rc3_data);
  tree_mc->SetBranchAddress("rc3",&rc3_mc);
  tree_data->SetBranchAddress("rc4",&rc4_data);
  tree_mc->SetBranchAddress("rc4",&rc4_mc);
  tree_mc->SetBranchAddress("weight",&weight);

  TH1 *h_rc_data = new TH1F("h_rc_data",";;",100,-10,10);
  TH1 *h_rc_mc = new TH1F("h_rc_mc",";;",100,-10,10);
  TH1 *h_rc_best_fit = new TH1F("h_rc_mc_best_fit",";;",100,-10,10);

  float rc_data=-1,rc_mc=-1;
  //data
  for (int e=0;e<tree_data->GetEntries();e++){
    tree_data->GetEntry(e);
    if(rc_number==2) rc_data = rc2_data;
    if(rc_number==3) rc_data = rc3_data;
    if(nring_data==rc_number || nring_data==rc_number+1) h_rc_data->Fill(rc_data);
  }

  //mc
  for (int e=0;e<tree_mc->GetEntries();e++){
    tree_mc->GetEntry(e);
    if(rc_number==2) rc_mc = rc2_mc;
    if(rc_number==3) rc_mc = rc3_mc;
    if(nring_mc==rc_number || nring_mc==rc_number+1) h_rc_mc->Fill(rc_mc,weight);
  }
  //h_rc_data->Draw("E0");
  //h_rc_mc->Draw("hist same");
  float norm = h_rc_data->Integral()/h_rc_mc->Integral();
  cout << "mc normalization=" << norm << endl;

  TCanvas *c = new TCanvas("c1","",800,600);
  h_rc_data->SetLineWidth(2);
  h_rc_data->Draw("E0");
  h_rc_mc->Scale(norm);
  h_rc_mc->SetLineWidth(2);
  h_rc_mc->SetLineColor(2);
  h_rc_mc->Draw("hist same");

  if(calc_chi2){//calculate chi2
    int nbin = h_rc_data->GetNbinsX();
    float min_chi2=9999.,min_scale=0.,min_shift=0.;
    for (int a=0;a<100;a++){
      for(int b=0;b<100;b++){
        TH1 *h_rc_mc_scale_shift = new TH1F(Form("h_rc_mc_scale%d_shift%d",a,b),";;",100,-10,10);
        h_rc_mc_scale_shift->SetLineColor(4);
        float scale = 0.5 + 0.01*a;
        float shift = 0.01*b -0.5;
        //float scale = 1;
        //float shift = 0;
        for (int e=0;e<tree_mc->GetEntries();e++){
          tree_mc->GetEntry(e);
          if(rc_number==2) rc_mc = rc2_mc;
          if(rc_number==3) rc_mc = rc3_mc;
          if(nring_mc==rc_number || nring_mc==rc_number+1) 
            h_rc_mc_scale_shift->Fill(rc_mc*scale+shift,weight*norm);
        }
        cout << "scale/shift=" << scale << "/" << shift << endl;
        float chi2=0.;
        for(int bb=0;bb<nbin;bb++){
          float event_data = h_rc_data->GetBinContent(bb+1);
          float event_mc = h_rc_mc_scale_shift->GetBinContent(bb+1);
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
  if(rc_number==2){
    if(period==1){
      final_scale = 1.01;
      final_shift = -0.19;
    }
    if(period==2){
      final_scale = 1.0;
      final_shift = -0.17;
    }
    if(period==3){
      final_scale = 1.0;
      final_shift = -0.02;
    }
    if(period==4){
      final_scale = 0.99;
      final_shift = 0.03;
    }
  }
  if(rc_number==3){
    if(period==1){
      final_scale = 1.0;
      final_shift = -0.34;
    }
    if(period==2){
      final_scale = 1.02;
      final_shift = -0.02;
    }
    if(period==3){
      final_scale = 1.05;
      final_shift = 0;
    }
    if(period==4){
      final_scale = 0.98;
      final_shift = -0.02;
    }
  }
  for (int e=0;e<tree_mc->GetEntries();e++){
    tree_mc->GetEntry(e);
    if(rc_number==2) rc_mc = rc2_mc;
    if(rc_number==3) rc_mc = rc3_mc;
    if(nring_mc==rc_number || nring_mc==rc_number+1) 
      h_rc_mc_best_fit->Fill(rc_mc*final_scale+final_shift,weight*norm);
  }
  h_rc_mc_best_fit->SetLineColor(4);
  h_rc_mc_best_fit->SetLineWidth(2);
  h_rc_mc_best_fit->Draw("same hist");
  if(!calc_chi2) c->SaveAs(Form("hist/compare_data_mc_rc%d_before_after_fit_subgev_multiring_sk%d.pdf",rc_number,period));

}

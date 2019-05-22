void syst_ring_calc(){
  //setup
  string input_name = "fcmc_rc";
  string mode_name = "p_mumumu";
  int period = 4;

  gROOT->SetStyle("Plain");
  gStyle->SetOptStat(0);

  TFile *input = TFile::Open(Form("../output/%s.sk%d.mode_%s_total_box_tree.root",input_name.c_str(),period,mode_name.c_str()));
  TTree *tree = (TTree*) input->Get("osc_tuple");
  TFile *output = new TFile(Form("output/output_syst_ring_input_%s_mode_%s_sk%d.root",input_name.c_str(),mode_name.c_str(),period),"recreate");

  int nring;
  float rc1,rc2,rc3,rc4,weight,total_mass,total_mom;
  tree->SetBranchAddress("nring",&nring);
  tree->SetBranchAddress("rc1",&rc1);
  tree->SetBranchAddress("rc2",&rc2);
  tree->SetBranchAddress("rc3",&rc3);
  tree->SetBranchAddress("rc4",&rc4);
  tree->SetBranchAddress("weight",&weight);
  tree->SetBranchAddress("total_mass",&total_mass);
  tree->SetBranchAddress("total_mom",&total_mom);

  TH1 *h_rc2_low = new TH1F("h_rc2_low",";;",100,-10,10);
  TH1 *h_rc2_high = new TH1F("h_rc2_high",";;",100,-10,10);
  TH1 *h_rc2_total = new TH1F("h_rc2_total",";;",100,-10,10);
  TH1 *h_rc3_low = new TH1F("h_rc3_low",";;",100,-10,10);
  TH1 *h_rc3_high = new TH1F("h_rc3_high",";;",100,-10,10);
  TH1 *h_rc3_total = new TH1F("h_rc3_total",";;",100,-10,10);
  TH1 *h_tuned_rc2_low = new TH1F("h_tuned_rc2_low",";;",100,-10,10);
  TH1 *h_tuned_rc2_high = new TH1F("h_tuned_rc2_high",";;",100,-10,10);
  TH1 *h_tuned_rc2_total = new TH1F("h_tuned_rc2_total",";;",100,-10,10);
  TH1 *h_tuned_rc3_low = new TH1F("h_tuned_rc3_low",";;",100,-10,10);
  TH1 *h_tuned_rc3_high = new TH1F("h_tuned_rc3_high",";;",100,-10,10);
  TH1 *h_tuned_rc3_total = new TH1F("h_tuned_rc3_total",";;",100,-10,10);

  float scale_rc2,shift_rc2,scale_rc3,shift_rc3;
  if(period==1){
    scale_rc2 = 1.01;
    shift_rc2 = -0.19;
    scale_rc3 = 1.0;
    shift_rc3 = -0.34;
  }
  if(period==2){
    scale_rc2 = 1.0;
    shift_rc2 = -0.17;
    scale_rc3 = 1.02;
    shift_rc3 = -0.02;
  }
  if(period==3){
    scale_rc2 = 1.0;
    shift_rc2 = -0.02;
    scale_rc3 = 1.05;
    shift_rc3 = 0;
  }
  if(period==4){
    scale_rc2 = 0.99;
    shift_rc2 = 0.03;
    scale_rc3 = 0.98;
    shift_rc3 = -0.02;
  }

  //mc
  float rc2_pos_total=0,tuned_rc2_pos_total=0,rc3_neg_total=0,tuned_rc3_neg_total=0;
  float rc2_pos_low=0,tuned_rc2_pos_low=0,rc3_neg_low=0,tuned_rc3_neg_low=0;
  float rc2_pos_high=0,tuned_rc2_pos_high=0,rc3_neg_high=0,tuned_rc3_neg_high=0;
  float err_rc2_pos_total=0,err_tuned_rc2_pos_total=0,err_rc3_neg_total=0,err_tuned_rc3_neg_total=0;
  float err_rc2_pos_low=0,err_tuned_rc2_pos_low=0,err_rc3_neg_low=0,err_tuned_rc3_neg_low=0;
  float err_rc2_pos_high=0,err_tuned_rc2_pos_high=0,err_rc3_neg_high=0,err_tuned_rc3_neg_high=0;
  float total_events=0;
  for (int e=0;e<tree->GetEntries();e++){
    tree->GetEntry(e);
    total_events+=weight;
    float tuned_rc2 = rc2*scale_rc2 + shift_rc2;
    float tuned_rc3 = rc3*scale_rc3 + shift_rc3;
    //cout << "rc1/rc2/rc3/rc4=" << rc1 << "/" << rc2 << "/" << rc3 << "/" << rc4 << endl;
    h_rc2_total->Fill(rc2,weight);
    h_rc3_total->Fill(rc3,weight);
    h_tuned_rc2_total->Fill(tuned_rc2,weight);
    h_tuned_rc3_total->Fill(tuned_rc3,weight);
    if(rc2>0) {
      rc2_pos_total+=weight;
      err_rc2_pos_total+=weight*weight;
    }
    if(rc3<0) {
      rc3_neg_total+=weight;
      err_rc3_neg_total+=weight*weight;
    }
    if(tuned_rc2>0) {
      tuned_rc2_pos_total+=weight;
      err_tuned_rc2_pos_total+=weight*weight;
    }
    if(tuned_rc3<0) {
      tuned_rc3_neg_total+=weight;
      err_tuned_rc3_neg_total+=weight*weight;
    }
    if(total_mom<100){//low signal box
      h_rc2_low->Fill(rc2,weight);
      h_rc3_low->Fill(rc3,weight);
      h_tuned_rc2_low->Fill(tuned_rc2,weight);
      h_tuned_rc3_low->Fill(tuned_rc3,weight);
      if(rc2>0) {
        rc2_pos_low+=weight;
        err_rc2_pos_low+=weight*weight;
      }
      if(rc3<0) {
        rc3_neg_low+=weight;
        err_rc3_neg_low+=weight*weight;
      }
      if(tuned_rc2>0) {
        tuned_rc2_pos_low+=weight;
        err_tuned_rc2_pos_low+=weight*weight;
      }
      if(tuned_rc3<0) {
        tuned_rc3_neg_low+=weight;
        err_tuned_rc3_neg_low+=weight*weight;
      }
    }
    else if(total_mom<250){//high signal box
      h_rc2_high->Fill(rc2,weight);
      h_rc3_high->Fill(rc3,weight);
      h_tuned_rc2_high->Fill(tuned_rc2,weight);
      h_tuned_rc3_high->Fill(tuned_rc3,weight);
      if(rc2>0) {
        rc2_pos_high+=weight;
        err_rc2_pos_high+=weight*weight;
      }
      if(rc3<0) {
        rc3_neg_high+=weight;
        err_rc3_neg_high+=weight*weight;
      }
      if(tuned_rc2>0) {
        tuned_rc2_pos_high+=weight;
        err_tuned_rc2_pos_high+=weight*weight;
      }
      if(tuned_rc3<0) {
        tuned_rc3_neg_high+=weight;
        err_tuned_rc3_neg_high+=weight*weight;
      }
    }
  }
  err_rc2_pos_total = sqrt(err_rc2_pos_total);
  err_tuned_rc2_pos_total = sqrt(err_tuned_rc2_pos_total);
  err_rc3_neg_total = sqrt(err_rc3_neg_total);
  err_tuned_rc3_neg_total = sqrt(err_tuned_rc3_neg_total);
  err_rc2_pos_low = sqrt(err_rc2_pos_low);
  err_tuned_rc2_pos_low = sqrt(err_tuned_rc2_pos_low);
  err_rc3_neg_low = sqrt(err_rc3_neg_low);
  err_tuned_rc3_neg_low = sqrt(err_tuned_rc3_neg_low);
  err_rc2_pos_high = sqrt(err_rc2_pos_high);
  err_tuned_rc2_pos_high = sqrt(err_tuned_rc2_pos_high);
  err_rc3_neg_high = sqrt(err_rc3_neg_high);
  err_tuned_rc3_neg_high = sqrt(err_tuned_rc3_neg_high);
  float diff_rc2_pos_total = tuned_rc2_pos_total - rc2_pos_total;
  float diff_rc3_neg_total = tuned_rc3_neg_total - rc3_neg_total;
  float ratio_rc2_pos_total = (rc2_pos_total>1e-5)? diff_rc2_pos_total / rc2_pos_total:0;
  float ratio_rc3_neg_total = (rc3_neg_total>1e-5)? diff_rc3_neg_total / rc3_neg_total:0;
  float stat_err_rc2_pos_total = (rc2_pos_total>1e-5)? err_rc2_pos_total/rc2_pos_total : 0;
  float stat_err_rc3_neg_total = (rc3_neg_total>1e-5)? err_rc3_neg_total/rc3_neg_total : 0;
  float stat_err_rc2_pos_low = (rc2_pos_low>1e-5)? err_rc2_pos_low/rc2_pos_low : 0;
  float stat_err_rc3_neg_low = (rc3_neg_low>1e-5)? err_rc3_neg_low/rc3_neg_low : 0;
  float stat_err_rc2_pos_high = (rc2_pos_high>1e-5)? err_rc2_pos_high/rc2_pos_high : 0;
  float stat_err_rc3_neg_high = (rc3_neg_high>1e-5)? err_rc3_neg_high/rc3_neg_high : 0;
  cout << "input:" << input_name << " mode:" << mode_name << " sk" << period << endl;
  cout << "Events in signal box " << total_events << endl;
  //cout << "total rc2>0 def/tune/ratio=" << rc2_pos_total << "+-" << err_rc2_pos_total << "/" << tuned_rc2_pos_total << "+-" << err_tuned_rc2_pos_total << "/" << ratio_rc2_pos_total << "+-" << ratio_rc2_pos_total*stat_err_rc2_pos_total << endl; 
  cout << "total rc3<0 def/tune/ratio=" << rc3_neg_total << "+-" << err_rc3_neg_total << "/" << tuned_rc3_neg_total << "+-" << err_tuned_rc3_neg_total << "/" << ratio_rc3_neg_total << "+-" << ratio_rc3_neg_total*stat_err_rc3_neg_total << endl; 
  float diff_rc2_pos_low = tuned_rc2_pos_low - rc2_pos_low;
  float diff_rc3_neg_low = tuned_rc3_neg_low - rc3_neg_low;
  float ratio_rc2_pos_low = (rc2_pos_low>1e-5)? diff_rc2_pos_low / rc2_pos_low:0;
  float ratio_rc3_neg_low = (rc3_neg_low>1e-5)? diff_rc3_neg_low / rc3_neg_low:0;
  //cout << "low rc2>0 def/tune/ratio=" << rc2_pos_low << "+-" << err_rc2_pos_low << "/" << tuned_rc2_pos_low << "+-" << err_tuned_rc2_pos_low << "/" << ratio_rc2_pos_low << "+-" << ratio_rc2_pos_low*stat_err_rc2_pos_low << endl; 
  cout << "low rc3<0 def/tune/ratio=" << rc3_neg_low << "+-" << err_rc3_neg_low << "/" << tuned_rc3_neg_low << "+-" << err_tuned_rc3_neg_low << "/" << ratio_rc3_neg_low << "+-" << ratio_rc3_neg_low*stat_err_rc3_neg_low << endl; 
  float diff_rc2_pos_high = tuned_rc2_pos_high - rc2_pos_high;
  float diff_rc3_neg_high = tuned_rc3_neg_high - rc3_neg_high;
  float ratio_rc2_pos_high = (rc2_pos_high>1e-5)? diff_rc2_pos_high / rc2_pos_high:0;
  float ratio_rc3_neg_high = (rc3_neg_high>1e-5)? diff_rc3_neg_high / rc3_neg_high:0;
  //cout << "high rc2>0 def/tune/ratio=" << rc2_pos_high << "+-" << err_rc2_pos_high << "/" << tuned_rc2_pos_high << "+-" << err_tuned_rc2_pos_high << "/" << ratio_rc2_pos_high << "+-" << ratio_rc2_pos_high*stat_err_rc2_pos_high << endl; 
  cout << "high rc3<0 def/tune/ratio=" << rc3_neg_high << "+-" << err_rc3_neg_high << "/" << tuned_rc3_neg_high << "+-" << err_tuned_rc3_neg_high << "/" << ratio_rc3_neg_high << "+-" << ratio_rc3_neg_high*stat_err_rc3_neg_high << endl; 

  output->Write();

}

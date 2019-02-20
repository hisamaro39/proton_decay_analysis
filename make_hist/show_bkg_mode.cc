void show_bkg_mode(){
  string type[] = {"p_epi","p_mupi","p_eee","p_mumumu","p_emumu","p_muee","p_eemu","fcmc","fcdt"};
  int pi0cut = 1;
  int mode_id = 5;
  int sk_period = 4;
  int cut_number = 7;
  //1:FC&FC 2:nRing 3:PID 4:decayE 5:pi0 mass 6:total mass&mom 7:ntag
  int nring=1;
  int mulike=0;
  int michel=1;
  cout << "mode: " << type[mode_id] << endl;
  cout << "cut type nring" << nring << " mulike" << mulike << " michel" << michel << endl;

  gROOT->SetStyle("Plain");
  gStyle->SetOptStat(0);
  gStyle->SetPalette(1);

  if(pi0cut) input = TFile::Open(Form("../output/%s.sk%d.mode_%s_pi0cut_check_bkg.root",type[7].c_str(),sk_period,type[mode_id].c_str()));
  else input = TFile::Open(Form("../output/%s.sk%d.mode_%s_check_bkg.root",type[7].c_str(),sk_period,type[mode_id].c_str()));
  TH1 *hist;
  int sort_type[2*2*55],sort_sign[2*2*55],sort_mode[2*2*55];
  float sort_events[2*2*55],sort_err[2*2*55];
  for(int i=0;i<2*2*55;i++){
    sort_type[i]=-1;
    sort_sign[i]=0;
    sort_mode[i]=0;
    sort_events[i]=0.;
    sort_err[i]=0.;
  }

  float total_events=0.,total_err_square=0.;
  for(int p=0;p<2;p++){
    for (int t=0;t<2;t++){
      for (int m=0;m<54;m++){
        if(p==0) hist = (TH1*) input->Get(Form("cut_flow_nring%d_mulike%d_michel%d_type%d_mode_pos%d",nring,mulike,michel,t,m+1)); 
        if(p==1) hist = (TH1*) input->Get(Form("cut_flow_nring%d_mulike%d_michel%d_type%d_mode_neg%d",nring,mulike,michel,t,m+1)); 
        float events = hist->GetBinContent(cut_number+1);
        float err = hist->GetBinError(cut_number+1);
        if(cut_number==6) {
          events = hist->GetBinContent(cut_number+1) + hist->GetBinContent(cut_number+2);
          err = sqrt( pow(hist->GetBinError(cut_number+1),2) + pow(hist->GetBinError(cut_number+2),2) );
        }
        if(cut_number==7) {
          events = hist->GetBinContent(cut_number+2) + hist->GetBinContent(cut_number+3);
          err = sqrt( pow(hist->GetBinError(cut_number+2),2) + pow(hist->GetBinError(cut_number+3),2) );
        }
        //cout << "events is " << events << " +- " << err << endl;
        total_events += events;
        total_err_square += pow(err,2);
        for(int s=0;s<2*2*55;s++){
          if(events > sort_events[s]){
            for(int r=2*2*55-1;r>s;r--){
              sort_events[r] = sort_events[r-1];
              sort_err[r] = sort_err[r-1];
              sort_type[r] = sort_type[r-1];
              sort_sign[r] = sort_sign[r-1];
              sort_mode[r] = sort_mode[r-1];
            }
            sort_events[s] = events;
            sort_err[s] = err;
            sort_type[s] = t;
            sort_sign[s] = p;
            sort_mode[s] = m+1;
            break;
          }
        }
      }
    }
  }
  float total_ratio=0.;
  TH1 *final_hist,*final_hist2;
  cout << "total_events=" << total_events << " +- " << sqrt(total_err_square) << endl;
  cout << "### sort result ###" << endl;
  for(int f=0;f<2*2*55;f++){
    if(sort_events[f]<1e-5) break;
    if(sort_type[f]==0) cout << "electron type ";
    else cout << "muon type ";
    int final_mode = (sort_sign[f]==0)? sort_mode[f] : -1*sort_mode[f];
    float ratio = 1.*sort_events[f]/total_events;
    if(sort_sign[f]==0) final_hist = (TH1*) input->Get(Form("mass_proton_reco_cut%d_nring%d_mulike%d_michel%d_type%d_mode_pos%d",cut_number,nring,mulike,michel,sort_type[f],sort_mode[f])); 
    else final_hist = (TH1*) input->Get(Form("mass_proton_reco_cut%d_nring%d_mulike%d_michel%d_type%d_mode_neg%d",cut_number,nring,mulike,michel,sort_type[f],sort_mode[f])); 
    int raw_events = final_hist->GetEntries();
    if(cut_number==6){
      if(sort_sign[f]==0) final_hist2 = (TH1*) input->Get(Form("mass_proton_reco_cut%d_nring%d_mulike%d_michel%d_type%d_mode_pos%d",cut_number+1,nring,mulike,michel,sort_type[f],sort_mode[f])); 
      else final_hist2 = (TH1*) input->Get(Form("mass_proton_reco_cut%d_nring%d_mulike%d_michel%d_type%d_mode_neg%d",cut_number+1,nring,mulike,michel,sort_type[f],sort_mode[f])); 
      raw_events += final_hist2->GetEntries();
    }
    if(cut_number==7){
      if(sort_sign[f]==0) {
        final_hist = (TH1*) input->Get(Form("mass_proton_reco_cut%d_nring%d_mulike%d_michel%d_type%d_mode_pos%d",cut_number+1,nring,mulike,michel,sort_type[f],sort_mode[f])); 
        final_hist2 = (TH1*) input->Get(Form("mass_proton_reco_cut%d_nring%d_mulike%d_michel%d_type%d_mode_pos%d",cut_number+2,nring,mulike,michel,sort_type[f],sort_mode[f])); 
      }
      else {
        final_hist = (TH1*) input->Get(Form("mass_proton_reco_cut%d_nring%d_mulike%d_michel%d_type%d_mode_neg%d",cut_number+1,nring,mulike,michel,sort_type[f],sort_mode[f])); 
        final_hist2 = (TH1*) input->Get(Form("mass_proton_reco_cut%d_nring%d_mulike%d_michel%d_type%d_mode_neg%d",cut_number+2,nring,mulike,michel,sort_type[f],sort_mode[f])); 
      }
      raw_events = final_hist->GetEntries() + final_hist2->GetEntries();
    }
    total_ratio+=ratio;
    cout << "mode: " << final_mode;
    cout << " events: " << sort_events[f] << " +- " << sort_err[f] << " raw: " << raw_events << " ratio: " << ratio << endl;
  }
  cout << "total ratio is " << total_ratio << endl;

}



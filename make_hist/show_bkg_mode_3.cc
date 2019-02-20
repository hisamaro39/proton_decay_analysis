void show_bkg_mode_3(){
  string type[] = {"subgev_multiring"};
  int mode_id = 0;
  int sk_period = 4;
  int n_michel = -1;
  int n_mulike = 1;
  cout << "mode: " << type[mode_id] << endl;
  cout << "check n_michel_electron=" << n_michel << endl;

  gROOT->SetStyle("Plain");
  gStyle->SetOptStat(0);
  gStyle->SetPalette(1);

  input = TFile::Open(Form("../output/fcmc.sk%d.mode_%s.root",sk_period,type[mode_id].c_str()));
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
        if(p==0) hist = (TH1*) input->Get(Form("n_michel_electron_mulike%d_type%d_mode_pos%d",n_mulike,t,m+1)); 
        if(p==1) hist = (TH1*) input->Get(Form("n_michel_electron_mulike%d_type%d_mode_neg%d",n_mulike,t,m+1)); 
        float events=0.; 
        if(n_michel==-1) events = hist->Integral();
        else events = hist->GetBinContent(n_michel+1);
        total_events += events;
        float err=0.;
        if(n_michel==-1) err = hist->Integral();
        else err = hist->GetBinError(n_michel+1);
        //cout << "events is " << events << " +- " << err << endl;
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
  float total_ratio=0.,sum_ratio=0.;
  TH1 *final_hist,*final_hist2;
  cout << "### sort result ###" << endl;
  for(int f=0;f<2*2*55;f++){
    //if(sort_events[f]<1e-5) break;
    if(sort_type[f]==0) cout << "electron type ";
    else cout << "muon type ";
    int final_mode = (sort_sign[f]==0)? sort_mode[f] : -1*sort_mode[f];
    float ratio = 1.*sort_events[f]/total_events;
    if(ratio<1e-2) break;
    sum_ratio += ratio;
    cout << "mode: " << final_mode;
    cout << " events: " << sort_events[f] << " +- " << sort_err[f] << " [" << ratio*100 << "%]" << endl;
  }
  cout << "sum ratio is " << sum_ratio << endl;

}



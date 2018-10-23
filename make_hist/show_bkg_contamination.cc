void show_bkg_contamination(){
  string type[] = {"p_epi","p_mupi","p_eee","p_mumumu","p_emumu","p_muee","fcmc","fcdt"};
  int mode_id = 2;
  int sk_period = 4;
  int cut_number = 7;
  int nring=1;
  int nring2=0;
  int mulike=0;
  int michel=0;
  cout << "mode: " << type[mode_id] << endl;
  cout << "cut type nring" << nring << " mulike" << mulike << " michel" << michel << endl;

  gROOT->SetStyle("Plain");
  gStyle->SetOptStat(0);
  gStyle->SetPalette(1);

  TFile *input = TFile::Open(Form("../output/%s.sk%d.mode_%s.root",type[6].c_str(),sk_period,type[mode_id].c_str()));
  TH1 *hist = (TH1*) input->Get(Form("interaction_mode_cut%d_nring%d_mulike%d_michel%d_fp0",cut_number,nring,mulike,michel)); 
  TH1 *hist2 = (TH1*) input->Get(Form("interaction_mode_cut%d_nring%d_mulike%d_michel%d_fp0",cut_number,nring2,mulike,michel)); 
  hist->Add(hist2);
  float total = hist->Integral();
  cout << "total events is " << total << endl;
  int sort_id[93];
  float sort_ratio[93];
  for(int i=0;i<93;i++){
    sort_id[i]=-1;
    sort_ratio[i]=0;
  }
  for(int b=0;b<93;b++){
    float ratio = hist->GetBinContent(b+1) / total;
    //cout << "ratio=" << ratio << endl;
    for(int s=0;s<93;s++){
      //cout << "sort_ratio=" << sort_ratio[s] << endl;
      if(ratio > sort_ratio[s]){
        //cout << "ratio > sort_ratio!!" << endl;
        for(int r=92;r>s;r--){
          sort_ratio[r] = sort_ratio[r-1];
          sort_id[r] = sort_id[r-1];
        }
        sort_ratio[s] = ratio;
        sort_id[s] = b+1;
        break;
      }
    }
  }
  for(int bb=0;bb<93;bb++) {
    if(sort_ratio[bb]<0.01) continue;
    cout << "nuType:" << sort_id[bb]/23 + 1 << " ID:" << sort_id[bb]%23 << " " << sort_ratio[bb]*100 << " %" << endl;
  }

}



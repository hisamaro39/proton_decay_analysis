#include <vector>
void show_cut_flow_combine(){
  string type[] = {"p_epi","p_mupi","p_eee","p_mumumu","p_emumu","p_muee","fcmc","fcdt"};
  //input & mode
  int input_type=6;
  int mode_type=5;
  int sk_period=4;
  //cut pattern
  vector<int> nring;//0:2ring 1:3ring 2:4ring
  vector<int> nmulike;//# of mulike ring
  vector<int> nmichel;//# of michel electron
  nring.clear();nmulike.clear();nmichel.clear();
  //# of ring
  //nring.push_back(0);
  nring.push_back(1);
  //# of mu-like ring
  nmulike.push_back(1);
  //# of decay electron
  nmichel.push_back(1);

  TFile *input = TFile::Open(Form("../output/%s.sk%d.mode_%s.root",type[input_type].c_str(),sk_period,type[mode_type].c_str()));
  const int ring_pattern = nring.size();
  const int mulike_pattern = nmulike.size();
  const int michel_pattern = nmichel.size();
  TH1 *cut_flow[ring_pattern][mulike_pattern][michel_pattern];
  cout << "#############################" << endl;
  cout << "Showing cut flow of input:" << type[input_type] << " mode:" << type[mode_type] << " sk" << sk_period << endl;
  cout << "Combine...." << endl;
  for(int r=0;r<nring.size();r++){ 
    for(int mu=0;mu<nmulike.size();mu++) {
      for(int m=0;m<nmichel.size();m++) {
        cout << "ring pattern:" << nring[r] << " mu-like pattern:" << nmulike[mu] << " michel pattern:" << nmichel[m] << endl; 
        cut_flow[r][mu][m] = (TH1*) input->Get(Form("cut_flow_nring%d_mulike%d_michel%d",nring[r],nmulike[mu],nmichel[m]));  
      }
    }
  }
  cout << "#############################" << endl;

  float first_yield = cut_flow[0][0][0]->GetBinContent(1);
  float previous_yield = first_yield;
  for(int b=0;b<9;b++){
    float this_yield=0.;

    if(b<2) this_yield = cut_flow[0][0][0]->GetBinContent(b+1);
    if(b==2){ 
      for(int r=0;r<ring_pattern;r++)
        this_yield += cut_flow[r][0][0]->GetBinContent(b+1);
    }
    if(b==3){ 
      for(int r=0;r<ring_pattern;r++)
        for(int mu=0;mu<mulike_pattern;mu++)
        this_yield += cut_flow[r][mu][0]->GetBinContent(b+1);
    }
    if(b>3){ 
      for(int r=0;r<ring_pattern;r++)
        for(int mu=0;mu<mulike_pattern;mu++)
          for(int m=0;m<nmichel.size();m++) 
            this_yield += cut_flow[r][mu][m]->GetBinContent(b+1);
    }

    float eff = this_yield / first_yield * 100.;
    float eff_to_prev = this_yield / previous_yield * 100.;
    cout << "cut" << b+1 << " eff. " << eff << "% ======= " 
      << " eff to prev. " << eff_to_prev << "% ======= "
      << " yield. " << this_yield << "evt" << endl;
    if(b==7) continue;
    previous_yield = this_yield;
  }

}

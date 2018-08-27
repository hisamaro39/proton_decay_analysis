#include <vector>
void show_cut_flow(){
  string type[] = {"p_epi","p_mupi","p_eee","p_mumumu","p_emumu","p_muee","fcmc","fcdt"};
  //input & mode
  int input_type=1;
  int mode_type=1;
  //cut pattern
  int nring=0;//0:2ring 1:3ring 2:4ring
  int nmulike=0;//# of mulike ring
  int nmichel=0;//# of michel electron

  TFile *input = TFile::Open(Form("../output/%s.sk4.mode_%s.root",type[input_type].c_str(),type[mode_type].c_str()));
  TH1* this_hist = (TH1*) input->Get(Form("cut_flow_nring%d_mulike%d_michel%d",nring,nmulike,nmichel));  
  int nbins = this_hist->GetNbinsX();
  float first_yield = this_hist->GetBinContent(1);
  float previous_yield = first_yield;
  cout << "Cut Flow input: " << type[input_type] << "   mode: " << type[mode_type] << endl; 
  cout << "Cut Pattern: nring" << nring << " nmulike" << nmulike << " nmichel" << nmichel << endl;
  for(int b=0;b<8;b++){
    float this_yield = this_hist->GetBinContent(b+1);
    float eff = this_yield / first_yield * 100.;
    float eff_to_prev = this_yield / previous_yield * 100.;
    cout << "cut" << b+1 << " eff. " << eff << "% ======= " 
      << " eff to prev. " << eff_to_prev << "% ======= "
      << " yield. " << this_yield << "evt" << endl;
    if(b==6) continue;
    previous_yield = this_yield;
  }
}

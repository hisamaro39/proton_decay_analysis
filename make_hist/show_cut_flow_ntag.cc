#include <vector>
void show_cut_flow_ntag(){
  string type[] = {"p_epi","p_mupi","p_eee","p_eee_miura","p_eee_final",//4
    "p_mumumu","p_mumumu_miura","p_mumumu_final","p_emumu","p_emumu_miura",//9
    "p_emumu_final","p_mumue_final","p_muee","p_muee_miura","p_muee_final",//14
    "p_eemu_final","p_eee","p_muee","fcmc_real","fcdt",//19
    "fcmc_final"};
  //input & mode
  int input_type=18;
  int mode_type=5;
  int period=4;//5:combine sk1-4
  int fp= (input_type<18)? 1 : 0;
  //cut pattern
  int nring=0;//# of ring
  int nmulike=0;//# of mulike ring
  int nmichel=3;//# of michel electron

  TH1 *this_hist,*this_hist_free;
  TFile *input;
  if(period==5) input = TFile::Open(Form("../output/%s.sk1_4.mode_%s_wo_livetime.root",type[input_type].c_str(),type[mode_type].c_str()));
  else input = TFile::Open(Form("../output/%s.sk%d.mode_%s.root",type[input_type].c_str(),period,type[mode_type].c_str()));
  this_hist = (TH1*) input->Get(Form("cut_flow_nring%d_mulike%d_michel%d",nring,nmulike,nmichel));  
  if(fp==1) this_hist_free = (TH1*) input->Get(Form("cut_flow_nring%d_mulike%d_michel%d_fp1",nring,nmulike,nmichel));  
  int nbins = this_hist->GetNbinsX();
  float first_yield = this_hist->GetBinContent(1);
  float first_yield_err = this_hist->GetBinError(1);
  float previous_yield = first_yield;
  cout << "Cut Flow sk" << period <<  " input: " << type[input_type] << "   mode: " << type[mode_type] << endl; 
  cout << "Cut Pattern: nring" << nring << " nmulike" << nmulike << " nmichel" << nmichel << endl;
  float total_eff = 0.,total_yield=0.,total_eff_free = 0.,total_yield_free=0.;
  float total_eff_err_square = 0.,total_eff_err_free_square=0.,total_yield_err_square=0.;
  float total_eff_nt = 0.,total_yield_nt=0.,total_eff_free_nt = 0.,total_yield_free_nt=0.;
  float total_eff_err_square_nt = 0.,total_eff_err_free_square_nt=0.,total_yield_err_square_nt=0.;
  for(int b=0;b<10;b++){
    float this_yield = this_hist->GetBinContent(b+1);
    float this_yield_err = this_hist->GetBinError(b+1);
    float this_yield_err_free = (fp==1)? this_hist_free->GetBinError(b+1) : 0;
    float eff = this_yield / first_yield * 100.;
    float eff_err = (first_yield)? sqrt(pow(this_yield*first_yield_err,2)+pow(first_yield*this_yield_err,2))/pow(first_yield,2) : 0;
    float this_yield_free = (fp==1)? this_hist_free->GetBinContent(b+1) : 0;
    float eff_free = this_yield_free / first_yield * 100.;
    float eff_err_free = (first_yield)? sqrt(pow(this_yield_free*first_yield_err,2)+pow(first_yield*this_yield_err_free,2))/pow(first_yield,2) : 0;
    float eff_to_prev = this_yield / previous_yield * 100.;
    if(b==6) cout << "Lower Signal Region without ntag" << endl;
    if(b==7) cout << "Higher Signal Region without ntag" << endl;
    if(b==8) cout << "Lower Signal Region with ntag" << endl;
    if(b==9) cout << "Higher Signal Region with ntag" << endl;
    cout << "cut" << b+1 << " eff. " << eff << " +- " << eff_err*100 << "% ======= " 
      << " eff to prev. " << eff_to_prev << "% ======= "
      << " yield. " << this_yield << " +- " << this_yield_err << "evt" << endl;
    cout << "cut" << b+1 << " free eff. " << eff_free << " +- " << eff_err_free*100 << "% ======= " 
      << " yield. " << this_yield_free << " +- " << this_yield_err_free << "evt" << endl;
    //if(b==7) continue;
    previous_yield = this_yield;
    if(b==6 || b==7){//without neutron tagging
      total_eff+=eff;
      total_yield+=this_yield;
      total_eff_free+=eff_free;
      total_yield_free+=this_yield_free;
      total_eff_err_square += pow(eff_err,2);
      total_yield_err_square += pow(this_yield_err,2);
    }
    if(b==8 || b==9){//with neutron tagging
      total_eff_nt+=eff;
      total_yield_nt+=this_yield;
      total_eff_free_nt+=eff_free;
      total_yield_free_nt+=this_yield_free;
      total_eff_err_square_nt += pow(eff_err,2);
      total_yield_err_square_nt += pow(this_yield_err,2);
    }
  }
  cout << "total without ntag eff. " << total_eff << " +- " << sqrt(total_eff_err_square)*100 << "% ======= " 
    << " yield. " << total_yield << " +- " << sqrt(total_yield_err_square) << "evt" << endl;
  cout << "total with ntag eff. " << total_eff_nt << " +- " << sqrt(total_eff_err_square_nt)*100 << "% ======= " 
    << " yield. " << total_yield_nt << " +- " << sqrt(total_yield_err_square_nt) << "evt" << endl;
  if(fp==1){
    cout << "free proton" << endl;
    cout << "total without ntag eff. " << total_eff_free << "% ======= " 
      << " yield. " << total_yield_free << "evt" << endl;
    cout << "total with ntag eff. " << total_eff_free_nt << "% ======= " 
      << " yield. " << total_yield_free_nt << "evt" << endl;
  }
}

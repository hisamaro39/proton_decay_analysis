#include <vector>
void show_cut_flow_sk4(){
  string type[] = {"p_epi","p_mupi","p_eee","p_mumumu","p_emumu","p_muee","p_eemu","fcmc","fcdt"};
  //input & mode
  int input_type=7;
  int mode_type=5;
  int period=4;
  int fp= (input_type<7)? 1 : 0;
  //cut pattern
  int nring=2;//0:2ring 1:3ring 2:4ring
  int nmulike=1;//# of mulike ring
  int nmichel=1;//# of michel electron

  TH1 *this_hist,*this_hist_free;
  TFile *input = TFile::Open(Form("../output/%s.sk%d.mode_%s.root",type[input_type].c_str(),period,type[mode_type].c_str()));
  this_hist = (TH1*) input->Get(Form("cut_flow_nring%d_mulike%d_michel%d",nring,nmulike,nmichel));  
  if(fp==1) this_hist_free = (TH1*) input->Get(Form("cut_flow_nring%d_mulike%d_michel%d_fp1",nring,nmulike,nmichel));  
  int nbins = this_hist->GetNbinsX();
  float first_yield = this_hist->GetBinContent(1);
  float first_yield_err = this_hist->GetBinError(1);
  float previous_yield = first_yield;
  cout << "Cut Flow input: " << type[input_type] << "   mode: " << type[mode_type] << endl; 
  cout << "Cut Pattern: nring" << nring << " nmulike" << nmulike << " nmichel" << nmichel << endl;
  float total_eff = 0.,total_yield=0.,total_eff_free = 0.,total_yield_free=0.;
  float total_eff_err_square = 0.,total_eff_err_free_square=0.,total_yield_err_square=0.;
  float total_eff_nt = 0.,total_yield_nt=0.,total_eff_free_nt = 0.,total_yield_free_nt=0.;
  float total_eff_err_square_nt = 0.,total_eff_err_free_square_nt=0.,total_yield_err_square_nt=0.;
  for(int b=0;b<10;b++){
    float this_yield = this_hist->GetBinContent(b+1);
    float this_yield_err = this_hist->GetBinError(b+1);
    float eff = this_yield / first_yield * 100.;
    float eff_err = (first_yield)? sqrt(pow(this_yield*first_yield_err,2)+pow(first_yield*this_yield_err,2))/pow(first_yield,2) : 0;
    float this_yield_free = (fp==1)? this_hist_free->GetBinContent(b+1) : 0;
    float eff_free = this_yield_free / first_yield * 100.;
    float eff_to_prev = this_yield / previous_yield * 100.;
    cout << "cut" << b+1 << " eff. " << eff << " +- " << eff_err*100 << "% ======= " 
      << " eff to prev. " << eff_to_prev << "% ======= "
      << " yield. " << this_yield << " +- " << this_yield_err << "evt" << endl;
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

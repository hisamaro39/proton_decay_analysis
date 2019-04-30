#include <vector>
void calc_syst(){
  //input & mode
  string input_type="fcmc_final";
  string mode_type="p_eee";
  string syst = "energydown_sk4";
  int period=5;//5:sk1-4
  int live=1;//0:w/o live time weight
  //cut pattern
  int nring=1;//# of ring
  int nmulike=0;//# of mulike ring
  int nmichel=0;//# of michel electron
  cout << "input is " << input_type << endl;
  cout << "mode is " << mode_type << endl;
  cout << "systematic is " << syst << endl;
  if(live==0) cout << "without live time weight" << endl;
  if(period==5) cout << "period is sk1-4" << endl; 
  else cout << "period is sk" << period << endl; 

  TH1 *this_hist,*this_hist_free;
  TFile *input,*input_syst;
  if(period==5){
    input = TFile::Open(Form("../output/%s.sk1_4.mode_%s_wo_livetime.root",input_type.c_str(),mode_type.c_str()));
    input_syst = TFile::Open(Form("../output/%s.sk1_4.mode_%s_%s.root",input_type.c_str(),mode_type.c_str(),syst.c_str()));
  }else {
    if(live==0) input = TFile::Open(Form("../output/%s.sk%d.mode_%s_wo_livetime.root",input_type.c_str(),period,mode_type.c_str()));
    else input = TFile::Open(Form("../output/%s.sk%d.mode_%s.root",input_type.c_str(),period,mode_type.c_str()));
    input_syst = TFile::Open(Form("../output/%s.sk%d.mode_%s_%s.root",input_type.c_str(),period,mode_type.c_str(),syst.c_str()));
  }
  TH1 *cut_flow = (TH1*) input->Get(Form("cut_flow_nring%d_mulike%d_michel%d",nring,nmulike,nmichel));  
  TH1 *cut_flow_syst = (TH1*) input_syst->Get(Form("cut_flow_nring%d_mulike%d_michel%d",nring,nmulike,nmichel));  
  float evt_low = cut_flow->GetBinContent(9);
  float evt_high = cut_flow->GetBinContent(10);
  float err_low = cut_flow->GetBinError(9);
  float err_high = cut_flow->GetBinError(10);
  float evt_total = evt_low + evt_high;
  float err_total = sqrt(pow(err_low,2)+pow(err_high,2));
  float evt_low_syst = cut_flow_syst->GetBinContent(9);
  float evt_high_syst = cut_flow_syst->GetBinContent(10);
  float err_low_syst = cut_flow_syst->GetBinError(9);
  float err_high_syst = cut_flow_syst->GetBinError(10);
  float evt_total_syst = evt_low_syst + evt_high_syst;
  float err_total_syst = sqrt(pow(err_low_syst,2)+pow(err_high_syst,2));
  cout << "# Default #" << endl;
  cout << "events in SR low/high=" << evt_low << "+-" << err_low << "/" << evt_high << "+-" << err_high << endl;
  cout << "events in total SR is " << evt_total << "+-" << err_total << endl; 
  cout << "# Systematic of " << syst << " #" << endl;
  cout << "events in SR low/high=" << evt_low_syst << "+-" << err_low_syst << "/" << evt_high_syst << "+-" << err_high_syst << endl;
  cout << "events in total SR is " << evt_total_syst << "+-" << err_total_syst << endl; 
  float ratio_low = (evt_low>1e-5)? (evt_low_syst-evt_low)/evt_low : 0;
  float stat_err_low = (evt_low>1e-5)? err_low/evt_low : 0;
  float ratio_high = (evt_high>1e-5)? (evt_high_syst-evt_high)/evt_high : 0;
  float stat_err_high = (evt_high>1e-5)? err_high/evt_high : 0;
  float ratio_total = (evt_total>1e-5)? (evt_total_syst-evt_total)/evt_total : 0;
  float stat_err_total = (evt_total>1e-5)? err_total/evt_total : 0;
  cout << "# Results #" << endl;
  cout << "Ratio low is " << ratio_low << "+-" << ratio_low*stat_err_low << endl;
  cout << "Ratio high is " << ratio_high << "+-" << ratio_high*stat_err_high << endl;
  cout << "Ratio total is " << ratio_total << "+-" << ratio_total*stat_err_total << endl;
}

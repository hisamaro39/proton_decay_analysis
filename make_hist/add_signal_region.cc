#include <vector>
void add_signal_region(){
  string type[] = {"p_epi","p_mupi","p_eee","p_mumumu","p_emumu","p_muee","fcmc","fcdt"};
  //input & mode
  int input_type=6;
  int mode_type=2;
  //cut pattern
  int nring=1;//# of ring
  int nmulike=0;//# of mulike ring
  int nmichel=0;//# of michel electron
  string syst = "default";
  //string syst = "energydown";
  //string syst = "energyup";
  //string syst = "nonunidown";
  //string syst = "nonuniup";

  cout << "systematic is " << syst << endl;
  TH1 *this_hist,*this_hist_free;
  TFile *input;
  int period=0;
  float total_events=0,total_err_squre=0;
  for(int p=0;p<4;p++){
    period = p+1;
    //if(p==0) period=1;
    //if(p==1) period=3;
    //if(p==2) period=4;
    cout << "SK-" << period << endl;
    if(syst=="default") input = TFile::Open(Form("../output/%s.sk%d.mode_%s.root",type[input_type].c_str(),period,type[mode_type].c_str()));
    else input = TFile::Open(Form("../output/%s.sk%d.mode_%s_%s.root",type[input_type].c_str(),period,type[mode_type].c_str(),syst.c_str()));
    this_hist = (TH1*) input->Get(Form("cut_flow_nring%d_mulike%d_michel%d",nring,nmulike,nmichel));  
    float event_sr_low = this_hist->GetBinContent(9);
    float event_sr_high = this_hist->GetBinContent(10);
    float err_sr_low = this_hist->GetBinError(9);
    float err_sr_high = this_hist->GetBinError(10);
    cout << "Events in SR low/high=" << event_sr_low << "/" << event_sr_high << endl;
    cout << "Events in total SR is " << event_sr_low + event_sr_high << " +- " 
      << sqrt(pow(err_sr_low,2)+pow(err_sr_high,2)) << endl;
    total_events += event_sr_low + event_sr_high;
    total_err_squre += pow(err_sr_low,2) + pow(err_sr_high,2);
  }
  float total_err = sqrt(total_err_squre);
  cout << "Event some of total period in SR is " << total_events << " +- " << total_err << endl;
}

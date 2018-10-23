#include <vector>
void show_significance(){
  string type[] = {"p_epi","p_mupi","p_eee","p_mumumu","p_emumu","p_muee"};
  //input & mode
  int mode_type=5;
  int sk_period=4;
  //cut pattern
  vector<int> nring;//0:2ring 1:3ring 2:4ring
  vector<int> nmulike;//# of mulike ring
  vector<int> nmichel;//# of michel electron

  TFile *input_signal = TFile::Open(Form("../output/%s.sk%d.mode_%s.root",type[mode_type].c_str(),sk_period,type[mode_type].c_str()));
  TFile *input_bkg = TFile::Open(Form("../output/fcmc.sk%d.mode_%s.root",sk_period,type[mode_type].c_str()));

  cout << "#############################" << endl;
  cout << "Showing signicicance of mode:" << type[mode_type] << " sk" << sk_period << endl;
  for(int ir=0;ir<3;ir++){ 
    for(int imu=0;imu<4;imu++) {
      for(int im=0;im<7;im++) {
        nring.clear();nmulike.clear();nmichel.clear();
        //# of ring
        if(ir==2) {
          nring.push_back(0);
          nring.push_back(1);
        }
        else nring.push_back(ir);
        //# of mu-like ring
        nmulike.push_back(imu);
        //# of decay electron
        if(im==4) {
          for(int te=0;te<2;te++) nmichel.push_back(te);
          //nmichel.push_back(0);
          //nmichel.push_back(1);
        }
        else if(im==5) for(int te=0;te<3;te++) nmichel.push_back(te);
        else if(im==6) for(int te=0;te<4;te++) nmichel.push_back(te);
        else nmichel.push_back(im);

        const int ring_pattern = nring.size();
        const int mulike_pattern = nmulike.size();
        const int michel_pattern = nmichel.size();
        TH1 *cut_flow_signal[4][4][4];
        TH1 *cut_flow_bkg[4][4][4];
        for(int r=0;r<ring_pattern;r++){ 
          for(int mu=0;mu<mulike_pattern;mu++) {
            for(int m=0;m<michel_pattern;m++) {
              cut_flow_bkg[r][mu][m] = (TH1*) input_bkg->Get(Form("cut_flow_nring%d_mulike%d_michel%d",nring[r],nmulike[mu],nmichel[m]));  
              cut_flow_signal[r][mu][m] = (TH1*) input_signal->Get(Form("cut_flow_nring%d_mulike%d_michel%d",nring[r],nmulike[mu],nmichel[m]));  
            }
          }
        }

        float first_yield_signal = cut_flow_signal[0][0][0]->GetBinContent(1);
        float final_signal_eff=0.,final_bkg_yield=0.;

        for(int b=0;b<9;b++){
          float yield_signal=0.,yield_bkg=0.;

          if(b<2) {
            yield_signal = cut_flow_signal[0][0][0]->GetBinContent(b+1);
            yield_bkg = cut_flow_bkg[0][0][0]->GetBinContent(b+1);
          }
          if(b==2){ 
            for(int r=0;r<ring_pattern;r++){
              yield_signal += cut_flow_signal[r][0][0]->GetBinContent(b+1);
              yield_bkg += cut_flow_bkg[r][0][0]->GetBinContent(b+1);
            }
          }
          if(b==3){ 
            for(int r=0;r<ring_pattern;r++){
              for(int mu=0;mu<mulike_pattern;mu++){
                yield_signal += cut_flow_signal[r][mu][0]->GetBinContent(b+1);
                yield_bkg += cut_flow_bkg[r][mu][0]->GetBinContent(b+1);
              }
            }
          }
          if(b>3){ 
            for(int r=0;r<ring_pattern;r++){
              for(int mu=0;mu<mulike_pattern;mu++){
                for(int m=0;m<nmichel.size();m++){ 
                  yield_signal += cut_flow_signal[r][mu][m]->GetBinContent(b+1);
                  yield_bkg += cut_flow_bkg[r][mu][m]->GetBinContent(b+1);
                }
              }
            }
          }

          float eff = (first_yield_signal > 1e-5)? yield_signal / first_yield_signal * 100. : 0;
          //cout << "cut" << b+1 << " signal eff. " << eff << "% ======= " 
            //<< " bkg yield. " << yield_bkg << "evt ====== " 
            //<< " significance is " << eff/sqrt(yield_bkg) << endl;

          if(b==7 || b==8){
            final_signal_eff += eff;
            final_bkg_yield += yield_bkg;
          }
        }
        if(final_signal_eff>0){
          float significance = final_signal_eff/sqrt(final_bkg_yield) ;
          float ses = (1 - final_bkg_yield)/(133.5*final_signal_eff);
          cout << "ring pattern:" << ir << " mu-like pattern:" << imu << " michel pattern:" << im << endl; 
          cout << "final signal_eff/bkg_yield/significance=" << final_signal_eff << "/" << final_bkg_yield << "/" << significance << endl;
          cout << "ses[(kton*years)^-1]=" << ses << endl;
        }
      }
    }
  }

}

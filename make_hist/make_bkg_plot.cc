#include <vector>
#include "THStack.h"
void make_bkg_plot(){
  string type[] = {"p_epi","p_mupi","p_eee","p_mumumu","p_emumu","p_muee","fcmc","fcdt"};

  gROOT->SetStyle("Plain");
  gStyle->SetOptStat(0);
  gStyle->SetPalette(1);
  vector<string> hist;
  vector<int> input_type,mode_type,sk_period;
  hist.clear();input_type.clear();mode_type.clear();sk_period.clear();
  //hist list
  hist.push_back("mass_proton_reco_cut6_nring0_mulike0_michel0");  
  input_type.push_back(6);mode_type.push_back(1);sk_period.push_back(4);

  TFile *input;
  for(int h=0;h<hist.size();h++){
    input = TFile::Open(Form("../output/%s.sk%d.mode_%s.root",type[input_type[h]].c_str(),sk_period[h],type[mode_type[h]].c_str()));
    //TCanvas *c = new TCanvas(Form("canvas%d",h),"",800,600);
    //string save_name = "";
    THStack *hs = new THStack("hs","");
    float total=0.;
    for(int m=1;m<93;m++){
      TH1* temp_hist = (TH1*) input->Get(Form("%s_mode%d",hist[h].c_str(),m));
      total += temp_hist->Integral();
    }
    float sum_ratio=0.;
    int color_type=2,other_id=0;
    TH1 *other_hist;
    for(int m=1;m<93;m++){
      TH1* this_hist = (TH1*) input->Get(Form("%s_mode%d",hist[h].c_str(),m));
      float ratio = this_hist->Integral()/total;
      if(ratio>0.05) {
        cout << "dominant mode!! " << m << endl;
        //cout << "ratio=" << ratio << endl;
        sum_ratio += ratio;
        this_hist->SetLineColor(color_type);
        this_hist->SetFillColor(color_type);
        this_hist->SetFillStyle(3001);
        hs->Add(this_hist);
        color_type++;
      }
      else{
        //cout << "other mode!! " << m << endl;
        //cout << "ratio=" << ratio << endl;
        if(other_id==0) {
          other_hist = this_hist;
          other_id=1;
        }
        else{
          other_hist->Add(this_hist);
        }
        //cout << "other ratio = " << other_hist->Integral()/total << endl;
      }
    }
    cout << "total=" << total << endl;
    cout << "sum_ratio=" << sum_ratio << endl;
    //cout << "other_ratio=" << other_hist->Integral()/total << endl;
    //other_hist->SetLineColor(color_type+1);
    //other_hist->SetFillColor(color_type+1);
    //other_hist->SetFillStyle(3001);
    //hs->Add(other_hist);
    hs->Draw("hist");
    //c->SaveAs(Form("hist/single_%s_input_%s_mode_%s.pdf",hist[h].c_str(),type[input_type[h]].c_str(),type[mode_type[h]].c_str()));
  }
}

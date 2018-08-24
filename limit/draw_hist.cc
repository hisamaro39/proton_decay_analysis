#include <vector>
void draw_hist(){

  gROOT->SetStyle("Plain");
  gStyle->SetOptStat(0);
  vector<string> hist;
  vector<vector<string>> hist_set;
  vector<int> scale,dology;
  hist.clear();hist_set.clear();scale.clear();dology.clear();

  TFile *input = TFile::Open("output_hist_limit/output.root");

  hist.push_back("bais_ncand0_eff0.2_syseff0.1_nbkg4_sysbkg0_factor63.7");
  hist.push_back("bais_ncand0_eff0.2_syseff0.1_nbkg8_sysbkg0_factor63.7");
  hist.push_back("bais_ncand0_eff0.2_syseff0.1_nbkg12_sysbkg0_factor63.7");
  scale.push_back(1);dology.push_back(0);
  hist_set.push_back(hist);hist.clear();  

  hist.push_back("bais_ncand0_eff0.2_syseff0.1_nbkg4_sysbkg0_factor63.7");
  hist.push_back("bais_ncand0_eff0.3_syseff0.1_nbkg4_sysbkg0_factor63.7");
  hist.push_back("bais_ncand0_eff0.4_syseff0.1_nbkg4_sysbkg0_factor63.7");
  scale.push_back(1);dology.push_back(0);
  hist_set.push_back(hist);hist.clear();  

  hist.push_back("bais_ncand0_eff0.2_syseff0.1_nbkg4_sysbkg0_factor63.7");
  hist.push_back("bais_ncand0_eff0.2_syseff0.1_nbkg4_sysbkg0.1_factor63.7");
  hist.push_back("bais_ncand0_eff0.2_syseff0.1_nbkg4_sysbkg0.2_factor63.7");
  scale.push_back(1);dology.push_back(0);
  hist_set.push_back(hist);hist.clear();  

  hist.push_back("bais_ncand0_eff0.2_syseff0.1_nbkg4_sysbkg0_factor63.7");
  hist.push_back("bais_ncand0_eff0.2_syseff0.2_nbkg4_sysbkg0_factor63.7");
  hist.push_back("bais_ncand0_eff0.2_syseff0.3_nbkg4_sysbkg0_factor63.7");
  scale.push_back(1);dology.push_back(0);
  hist_set.push_back(hist);hist.clear();  

  hist.push_back("bais_ncand1_eff0.2_syseff0.1_nbkg4_sysbkg0_factor63.7");
  hist.push_back("bais_ncand1_eff0.3_syseff0.1_nbkg4_sysbkg0_factor63.7");
  hist.push_back("bais_ncand1_eff0.4_syseff0.1_nbkg4_sysbkg0_factor63.7");
  scale.push_back(1);dology.push_back(0);
  hist_set.push_back(hist);hist.clear();  

  hist.push_back("bais_ncand1_eff0.2_syseff0.1_nbkg4_sysbkg0_factor63.7");
  hist.push_back("bais_ncand1_eff0.2_syseff0.1_nbkg8_sysbkg0_factor63.7");
  hist.push_back("bais_ncand1_eff0.2_syseff0.1_nbkg12_sysbkg0_factor63.7");
  scale.push_back(1);dology.push_back(0);
  hist_set.push_back(hist);hist.clear();

  hist.push_back("bais_ncand0_eff0.2_syseff0.1_nbkg4_sysbkg0_factor63.7");
  hist.push_back("bais_ncand1_eff0.2_syseff0.1_nbkg4_sysbkg0_factor63.7");
  scale.push_back(1);dology.push_back(0);
  hist_set.push_back(hist);hist.clear();  

  TFile *input;
  TH1 *first_hist;
  for(int s=0;s<hist_set.size();s++){
    float evtmax=0,evtmax_scale=0,evt_min=99999,evtmin_scale=99999;
    for(int ss=0;ss<hist_set[s].size();ss++){//decide max event of the hist
      TH1 *temp_hist = (TH1*) input->Get(hist_set[s].at(ss).c_str());
     // cout << "maximum=" << temp_hist->GetMaximum()/temp_hist->Integral() << endl;
      float maximum=0;
      for(int b=0;b<temp_hist->GetNbinsX();b++){
        //cout << "b/value=" << b << "/" << temp_hist->GetBinContent(b+1)/temp_hist->Integral() << endl;
        float this_value = temp_hist->GetBinContent(b+1)/temp_hist->Integral();
        if(this_value>maximum) maximum = this_value;
      }
      //cout << "maximum=" << maximum << endl;
      if(temp_hist->GetMaximum() > evtmax) evtmax = temp_hist->GetMaximum();
      if(maximum > evtmax_scale){ 
        evtmax_scale = maximum;
        evtmin_scale = 1./temp_hist->Integral();
      }
      if(ss==0) first_hist=temp_hist;
    }
    //cout << "evtmax_scale=" << evtmax_scale << endl;
    //cout << "evtmax=" << evtmax << endl;
    TCanvas *c = new TCanvas(Form("canvas%d",s),"",800,600);
    if(dology[s]) c->SetLogy();
    string save_name = "";
    save_name = "output_hist_limit/compare_";
    for(int h=0;h<hist_set[s].size();h++){
      TH1* this_hist = (TH1*) input->Get(hist_set[s].at(h).c_str());
      //calc 90% point
      float inte=0., rate_limit=0.;
      for(int b=0;b<this_hist->GetNbinsX();b++){
        inte += this_hist->GetBinContent(b+1);
        ratio = inte/this_hist->Integral();
        //cout << "ratio=" << ratio << endl;
        if(ratio>0.9){
          rate_limit = 0.002*(b+1);
          break;
        }
      }
      cout << "rate_limit=" << rate_limit << endl;
      if(scale[s] && this_hist->GetEntries()) this_hist->Scale(1./this_hist->Integral());
      this_hist->SetLineWidth(2);
      if(h==0){
        this_hist->GetXaxis()->SetRangeUser(0,0.3);
        if(scale[s]) this_hist->GetYaxis()->SetRangeUser(0,evtmax_scale*1.2);
        else this_hist->GetYaxis()->SetRangeUser(0,evtmax*1.2);
        if(dology[s]){
          if(scale[s]) this_hist->GetYaxis()->SetRangeUser(0.00001,evtmax_scale*1.2);
          else this_hist->GetYaxis()->SetRangeUser(0.00001,evtmax*1.2);
        }
        this_hist->SetLineColor(1);
        this_hist->Draw("hist");
        first_hist = this_hist;
      }
      else{
        this_hist->SetLineColor(h+1);
        this_hist->Draw("same hist");
      }
      TLine *line = new TLine(rate_limit,0,rate_limit,evtmax_scale*1.2);
      line->SetLineWidth(2);
      line->SetLineStyle(2);
      line->SetLineColor(h+1);
      line->Draw();
      save_name += hist_set[s].at(h);
    }
    save_name += ".pdf";
    c->SaveAs(save_name.c_str());
    //c->SaveAs("hist/temp.pdf");

  }


}

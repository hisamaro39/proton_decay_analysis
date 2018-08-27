#include <vector>
void combine_hist(){

  gROOT->SetStyle("Plain");
  gStyle->SetOptStat(0);
  vector<string> hist;
  vector<vector<string>> hist_set;
  hist.clear();hist_set.clear();

  TFile *input = TFile::Open("output_hist_limit/output.root");

  hist.push_back("bais_ncand0_eff0.1_syseff0.1_nbkg4_sysbkg0_factor63.7");
  hist.push_back("bais_ncand0_eff0.1_syseff0.1_nbkg4_sysbkg0_factor63.7");
  hist_set.push_back(hist);hist.clear();  

  TH1 *combined_hist; 
  for(int s=0;s<hist_set.size();s++){
    for(int h=0;h<hist_set[s].size();h++){
      TH1* this_hist = (TH1*) input->Get(hist_set[s].at(h).c_str());
      if(h==0) {
        combined_hist=this_hist;
        continue;
      }
      else {
        cout << "bin size=" << this_hist->GetNbinsX() << endl;
        for(int b=0;b<this_hist->GetNbinsX();b++){
          float new_value = combined_hist->GetBinContent(b+1) * this_hist->GetBinContent(b+1);
          combined_hist->SetBinContent(b+1, new_value);
        }
      }
    }
  }

  combined_hist->Scale(1./combined_hist->Integral());
  combined_hist->Draw();

}

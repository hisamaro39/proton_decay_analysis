#include <iostream>
#include <algorithm>
#include <unistd.h>

void make_hist(){
  int nbin=500;

  vector<string> table_names;
  table_names.clear();
  DIR* dp=opendir("./output_table");
  if (dp!=NULL)
  {
    struct dirent* dent;
    do{
      dent = readdir(dp);
      if (dent!=NULL) {
        char *str_tmp = dent->d_name;
        if(strstr(str_tmp,"txt")!=NULL){
          table_names.push_back(str_tmp);
        }
      }
    } while(dent!=NULL);
    closedir(dp);
  }

  cout << "size of file is " << table_names.size() << endl;

  TFile *output = new TFile("output_hist_limit/output.root","recreate");
  for(int t=0;t<table_names.size();t++){

    std::ifstream table(Form("./output_table/%s",table_names[t].c_str()));
    int length = table_names[t].length();
    string tmp = table_names[t].substr(0,length-4); 
    cout << "tmp=" << tmp << endl;
    TH1 *hist = new TH1F(tmp.c_str(),"",nbin,0,1);

    std::string line;
    double xbinlow,xbinhigh,value;
    for(int l=0;l<nbin;l++){//energy
      std::getline(table, line);
      sscanf(line.c_str(),"%lf %lf %lf",&xbinlow, &xbinhigh, &value);
      //cout << "xbinlow/xbinhigh/value=" << xbinlow << "/" << xbinhigh << "/" << value << endl;
      double mean = (xbinlow+xbinhigh)*0.5;
      hist->Fill(mean,value);
    }

  }
  output->Write();

}

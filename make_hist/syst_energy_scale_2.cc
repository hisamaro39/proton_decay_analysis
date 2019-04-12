#include <vector>
void syst_energy_scale_2(){
  string type[] = {"p_epi","p_mupi","p_eee_miura","p_muee","p_mumumu","p_emumu","p_mumue","p_muee","fcmc","fcdt"};
  int input_type = 2;
  int mode_type = 2;
  int cut = 5;
  int nring = 1;
  int mulike = 0;
  int michel = 0;

  gROOT->SetStyle("Plain");
  gStyle->SetOptStat(0);
  gStyle->SetPalette(1);
  //vector<string> graph;
  //vector<int> input_type,mode_type;
  //graph.clear();input_type.clear();mode_type.clear();
  //hist list
  //input_type.push_back(8);mode_type.push_back(2);

  TFile *input;
  input = TFile::Open(Form("../output/%s.sk4.mode_%s.root",type[input_type].c_str(),type[mode_type].c_str()));
  TCanvas *c = new TCanvas("canvas","",800,600);
  TH1* frame=gPad->DrawFrame(0, 0, 1250, 1000);
  string save_name = "";
  TGraph* graph = (TGraph*) input->Get(Form("mass_mom_proton_reco_ntag_cut%d_nring%d_mulike%d_michel%d",cut,nring,mulike,michel));  
  graph->Draw("p");
  //double *x = = graph->GetX();
  //cout << "size is " << x->size() << endl;
  int npoints = graph->GetN();
  cout << "# of point is " << npoints << endl;
  double mass,mom;
  int n_def_low=0,n_mass_up_low=0,n_mass_down_low=0,n_mom_up_low=0,n_mom_down_low=0;
  int n_def_high=0,n_mass_up_high=0,n_mass_down_high=0,n_mom_up_high=0,n_mom_down_high=0;
  for(int i=0;i<npoints;i++){
    graph->GetPoint(i,mass,mom);
    if(mass>800 && mass<1050 && mom<100) n_def_low++;
    if(mass>800 && mass<1050 && mom>100 && mom<250) n_def_high++;
    if(mass>800*1.021 && mass<1050*1.021 && mom<100) n_mass_up_low++;
    if(mass>800*1.021 && mass<1050*1.021 && mom>100 && mom<250) n_mass_up_high++;
    if(mass>800*0.979 && mass<1050*0.979 && mom<100) n_mass_down_low++;
    if(mass>800*0.979 && mass<1050*0.979 && mom>100 && mom<250) n_mass_down_high++;
    if(mass>800 && mass<1050 && mom<100*1.021) n_mom_up_low++;
    if(mass>800 && mass<1050 && mom>100*1.021 && mom<250*1.021) n_mom_up_high++;
    if(mass>800 && mass<1050 && mom<100*0.979) n_mom_down_low++;
    if(mass>800 && mass<1050 && mom>100*0.979 && mom<250*0.979) n_mom_down_high++;
  }
  float diff_mass_up_low = 1.*(n_mass_up_low - n_def_low)/n_def_low;
  float diff_mass_down_low = 1.*(n_mass_down_low - n_def_low)/n_def_low;
  float diff_mass_up_high = 1.*(n_mass_up_high - n_def_high)/n_def_high;
  float diff_mass_down_high = 1.*(n_mass_down_high - n_def_high)/n_def_high;
  float diff_mom_up_low = 1.*(n_mom_up_low - n_def_low)/n_def_low;
  float diff_mom_down_low = 1.*(n_mom_down_low - n_def_low)/n_def_low;
  float diff_mom_up_high = 1.*(n_mom_up_high - n_def_high)/n_def_high;
  float diff_mom_down_high = 1.*(n_mom_down_high - n_def_high)/n_def_high;
  cout << "n_def low/high=" << n_def_low << "/" << n_def_high << endl;
  cout << "n_mass_up low/high=" << n_mass_up_low << "/" << n_mass_up_high << endl;
  cout << "diff_mass_up low/high=" << diff_mass_up_low << "/" << diff_mass_up_high << endl;
  cout << "n_mass_down low/high=" << n_mass_down_low << "/" << n_mass_down_high << endl;
  cout << "diff_mass_down low/high=" << diff_mass_down_low << "/" << diff_mass_down_high << endl;
  cout << "n_mom_up low/high=" << n_mom_up_low << "/" << n_mom_up_high << endl;
  cout << "diff_mom_up low/high=" << diff_mom_up_low << "/" << diff_mom_up_high << endl;
  cout << "n_mom_down low/high=" << n_mom_down_low << "/" << n_mom_down_high << endl;
  cout << "diff_mom_down low/high=" << diff_mom_down_low << "/" << diff_mom_down_high << endl;

  /*TGraph* graph_sr_low = (TGraph*) input->Get(Form("mass_mom_proton_reco_cut6_nring%d_mulike%d_michel%d",nring,mulike,michel));  
  TGraph* graph_sr_high = (TGraph*) input->Get(Form("mass_mom_proton_reco_cut7_nring%d_mulike%d_michel%d",nring,mulike,michel));  
  TBox *box = new TBox(800,0,1050,100);
  box->SetFillStyle(0);
  box->SetLineWidth(2);
  box->SetLineColor(1);
  box->Draw();
  TBox *box2 = new TBox(800,100,1050,250);
  box2->SetFillStyle(0);
  box2->SetLineWidth(2);
  box2->SetLineColor(1);
  box2->Draw();
  TEllipse *e = new TEllipse(925, 125, 125, 125);
  e->SetFillStyle(0);
  e->SetLineColor(2);
  e->SetLineWidth(2);
  e->Draw();
  TEllipse *e2 = new TEllipse(925, 125, 500, 500);
  e2->SetFillStyle(0);
  e2->SetLineColor(2);
  e2->SetLineWidth(2);
  e2->SetLineStyle(2);
  e2->Draw();
  TEllipse *e3 = new TEllipse(925, 125, 1000, 1000);
  e3->SetFillStyle(0);
  e3->SetLineColor(2);
  e3->SetLineWidth(2);
  e3->SetLineStyle(2);
  e3->Draw();
  graph_all->Draw("p");
  graph_sr_low->SetMarkerStyle(8);
  graph_sr_low->Draw("p");
  graph_sr_high->SetMarkerStyle(8);
  graph_sr_high->Draw("p");
  c->SaveAs(Form("hist/single_mass_mom_proton_input_%s_mode_%s.pdf",type[input_type].c_str(),type[mode_type].c_str()));
  */

}

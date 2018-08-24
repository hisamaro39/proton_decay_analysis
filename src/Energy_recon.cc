#include <iostream>
#include <stdlib.h>
#include <math.h>
#include <string>

#include "TF1.h"
#include "TGraphErrors.h"
#include "TTree.h"
#include "TFile.h"
#include "TH2F.h"

#include "Energy_recon.h"
#include "osc_types.h"

const Int_t nbins_amom=55;
const Float_t amom_binwidth=0.1;

Float_t bin_value_median[nbins_type][nbins_nn][nbins_muedk];
Float_t bin_RMS[nbins_type][nbins_nn][nbins_muedk];

Energy_recon::Energy_recon(){
  os2 = new outputOscStructure;
}
Energy_recon::~Energy_recon(){
  delete os2;
}

void Energy_recon::SetOutputFile(TFile* file)
{
  ofile=file;
}

void Energy_recon::GetFunction(std::string func_name){
  TFile* file_in=new TFile(func_name.c_str());

  TTree* func_params=(TTree*)file_in->Get("func_params");
  Int_t bin_nn;
  Int_t bin_muedk;
  Int_t rec_type;
  Float_t param_value[5];
  Float_t offset_value;
  func_params->SetBranchAddress("rec_type",&rec_type);
  func_params->SetBranchAddress("bin_nn",&bin_nn);
  func_params->SetBranchAddress("bin_muedk",&bin_muedk);
  func_params->SetBranchAddress("param_value",param_value);
  func_params->SetBranchAddress("offset_value",&offset_value);
  Int_t nparms=func_params->GetEntries();
  for(Int_t i=0;i<nparms;i++){
    func_params->GetEntry(i);
    for(Int_t l=0;l<5;l++){
      fit_params[rec_type][bin_nn][bin_muedk][l]=param_value[l];
    }
    amom_offset[rec_type][bin_nn][bin_muedk]=offset_value;
  }
  file_in->Close();  
}

void Energy_recon::SaveFunction(std::string func_name){
  TFile* file_out=new TFile(func_name.c_str(),"recreate");
  TTree* func_params=new TTree("func_params","func_params");
  Int_t bin_nn;
  Int_t bin_muedk;
  Int_t rec_type;
  Float_t param_value[5];
  Float_t offset_value;
  func_params->Branch("rec_type",&rec_type,"rec_type/I");
  func_params->Branch("bin_nn",&bin_nn,"bin_nn/I");
  func_params->Branch("bin_muedk",&bin_muedk,"bin_muedk/I");
  func_params->Branch("param_value",param_value,"param_value[5]/F");
  func_params->Branch("offset_value",&offset_value,"offset_value/F");

  for(Int_t i=0;i<nbins_rec_type;i++){
    for(Int_t j=0;j<nbins_nn;j++){
      for(Int_t k=0;k<nbins_muedk;k++){
	rec_type=i;
	bin_nn=j;
	bin_muedk=k;
	for(Int_t l=0;l<5;l++){
	  param_value[l]=fit_params[rec_type][bin_nn][bin_muedk][l];
	}
	offset_value=amom_offset[rec_type][bin_nn][bin_muedk];
	func_params->Fill();
      }
    }
  }
  func_params->Write();
  file_out->Close();
}

void Energy_recon::GetSeparation(std::string cardname){
  ifstream ifs(cardname.c_str());
  int rec_type_id;
  int type;
  while(ifs>>type>>rec_type_id){
    rec_type[type-1]=rec_type_id;
  }
  for(Int_t i=0;i<nbins_type;i++){
    if(rec_type[i]<0||rec_type[i]>=nbins_rec_type) rec_type[i]=-1;
  }
    
  ifs.close();
}

void Energy_recon::ReadTree(OscNtupleManager* om){
  osc_tuple_origin=(TTree*)om->SendTree("osc_tuple_origin");
}

void Energy_recon::SendTree(OscNtupleManager* om){
  om->ReadTree(otree_erecon);
}

/*
void Energy_recon::DeleteTree(){
  otree_erecon->Delete();
}
*/
void Energy_recon::make_table(){

  TH2F* bin_value[nbins_rec_type][nbins_nn][nbins_muedk];
  
  Int_t nevent;
  Float_t amom;
  Float_t amom_pri;
  Float_t mass;
  Float_t pnu;
  Int_t nn=0;
  Int_t nn_mctruth=0;
  Int_t itype;
  Int_t muedk;
  Int_t ipnu;
  Float_t oscweight3f;
  Float_t dir[3];
  Float_t dirnu[3];
  Int_t mode;
  Float_t weightx;
  Float_t amom_forcorr;
  osc_tuple_origin->SetBranchAddress("amom",&amom);
  osc_tuple_origin->SetBranchAddress("amom_forcorr",&amom_forcorr);
  osc_tuple_origin->SetBranchAddress("amom_pri",&amom_pri);
  osc_tuple_origin->SetBranchAddress("pnu",&pnu);
  osc_tuple_origin->SetBranchAddress("nn",&nn);
  osc_tuple_origin->SetBranchAddress("nn_mctruth",&nn_mctruth);
  osc_tuple_origin->SetBranchAddress("itype",&itype);
  osc_tuple_origin->SetBranchAddress("muedk",&muedk);
  osc_tuple_origin->SetBranchAddress("ipnu",&ipnu);
  osc_tuple_origin->SetBranchAddress("oscweight3f",&oscweight3f);
  osc_tuple_origin->SetBranchAddress("weightx",&weightx);
  osc_tuple_origin->SetBranchAddress("dir",dir);
  osc_tuple_origin->SetBranchAddress("dirnu",dirnu);
  osc_tuple_origin->SetBranchAddress("mode",&mode);

  nevent=osc_tuple_origin->GetEntries();
  for(int bin_rec_type=0;bin_rec_type<nbins_rec_type;bin_rec_type++){    
    for(int bin_nn=0;bin_nn<nbins_nn;bin_nn++){    
      for(int bin_muedk=0;bin_muedk<nbins_muedk;bin_muedk++){
	bin_value[bin_rec_type][bin_nn][bin_muedk]=new TH2F(Form("bin_value_%d_%d_%d",bin_rec_type+1,bin_nn,bin_muedk),"",nbins_amom,1,1+nbins_amom*amom_binwidth,20000,-100000,100000);
      }
    }
  }
  Double_t xq[1]={0.5};
  Double_t yq[1];
  Double_t x[nbins_amom];	   
  Double_t y[nbins_amom];
  Int_t n=0;
  Double_t yerror[nbins_amom];
  TGraphErrors* g1 = new TGraphErrors(n,x,y,0,yerror);
  TF1* f1 = new TF1("f1","[0]+[1]*(pow(10,x)-pow(10,[3]))+([2]*(x>[3])+[4]*(x<[3]))*pow((pow(10,x)-pow(10,[3])),2)",0,6);
  for(Int_t i=0;i<5;i++){
    f1->SetParameter(i,0);
  }
  g1->Fit("f1","","",3.1,6);
  g1->Fit("f1","","",2.8,6);

  Int_t nnbin_input;
  Int_t muedkbin_input;
  Int_t rec_type_bin;
  Int_t bin_true[nbins_type][nbins_nn][nbins_muedk];
  for(Int_t i=0;i<nevent;i++){
    osc_tuple_origin->GetEntry(i);
    rec_type_bin=rec_type[(itype-1)%(EndOfTypes-1)];
    if(rec_type_bin==0||rec_type_bin==2) mass=0.511;    // Single-ring and Multi-ring elike samples
    else if(rec_type_bin==1||rec_type_bin==3) mass=105.6;    // Single-ring and Multi-ring mulike samples
    else mass=0;
    if(abs(mode)<=27){ //Use Only CC events
      if(nn<=6)nnbin_input=nn;
      else nnbin_input=(Int_t)6*(log((Float_t)nn/6.)/2.+1);

      if(muedk<=3)muedkbin_input=muedk;
      else muedkbin_input=muedk;

      if(rec_type_bin<0&&rec_type_bin>=nbins_rec_type) continue;
      if(nnbin_input>=nbins_nn) continue;
      if(muedkbin_input>=nbins_muedk) continue;
      bin_true[rec_type_bin][nnbin_input][muedkbin_input]++;
      if(rec_type_bin==0||rec_type_bin==1){ // Single Ring events
	bin_value[rec_type_bin][nnbin_input][muedkbin_input]->Fill(log10(amom_pri),(pnu*1000-sqrt(amom_pri*amom_pri+mass*mass)),weightx);
      }
      else if(rec_type_bin==2||rec_type_bin==3){ // Multi Ring events
	bin_value[rec_type_bin][nnbin_input][muedkbin_input]->Fill(log10(amom_pri),(pnu*1000-sqrt(amom_pri*amom_pri+mass*mass)-(amom_forcorr-sqrt(amom_pri*amom_pri+mass*mass))),weightx);
      }
    }  
  }

  TGraphErrors* g2[nbins_rec_type][nbins_nn][nbins_muedk];
  for(Int_t bin_rec_type=0;bin_rec_type<nbins_rec_type;bin_rec_type++){
    for(Int_t bin_nn=0;bin_nn<nbins_nn;bin_nn++){
      for(Int_t bin_muedk=0;bin_muedk<nbins_muedk;bin_muedk++){
	n=0;
	for(Int_t bin_amom=0;bin_amom<nbins_amom;bin_amom++){
	  TH1F* bin_value_py=(TH1F*)bin_value[bin_rec_type][bin_nn][bin_muedk]->ProjectionY("bin_value_py",bin_amom+1,bin_amom+1);
	  if(bin_value_py->Integral()>=5){
	    bin_value_py->GetQuantiles(1,yq,xq);
	    x[n]=bin_value[bin_rec_type][bin_nn][bin_muedk]->GetXaxis()->GetBinLowEdge(bin_amom+1)+bin_value[bin_rec_type][bin_nn][bin_muedk]->GetXaxis()->GetBinWidth(bin_amom+1)*0.5;
	  y[n]=yq[0];
	  yerror[n]=MAD(bin_value_py,yq[0])/sqrt((Double_t)bin_value_py->Integral());
	  bin_value_median[bin_rec_type][bin_nn][bin_muedk]=yq[0];
	  n++;
	}
	else{
	  bin_value_median[bin_rec_type][bin_nn][bin_muedk]=0;
	  bin_RMS[bin_rec_type][bin_nn][bin_muedk]=1e10;
	}
      }
      g2[bin_rec_type][bin_nn][bin_muedk]=new TGraphErrors(n,x,y,0,yerror);
      for(Int_t i=0;i<5;i++){
	f1->SetParameter(i,0);
      }
      if(bin_rec_type==3){ // MultiRing mu-like
	f1->SetParLimits(2,-1e-4,10);
	f1->FixParameter(3,3.1);
	f1->FixParameter(4,0);
	g2[bin_rec_type][bin_nn][bin_muedk]->Fit("f1","","",2.7,3.9);
	f1->SetParLimits(3,-1,10);
	g2[bin_rec_type][bin_nn][bin_muedk]->Fit("f1","","",2.7,3.9);
	g2[bin_rec_type][bin_nn][bin_muedk]->Fit("f1","","",2.7,3.9);
      }
      else{
	f1->FixParameter(2,0);
	f1->FixParameter(3,0);
	if(bin_rec_type==1){ // SingleRing mu-like
	  g2[bin_rec_type][bin_nn][bin_muedk]->Fit("f1","","",0,3.9);
	}
	else{
	  g2[bin_rec_type][bin_nn][bin_muedk]->Fit("f1","","",0,6);
	}	  
      }
      for(Int_t k=0;k<5;k++){
	if(n>=4)                  fit_params[bin_rec_type][bin_nn][bin_muedk][k]=f1->GetParameter(k);
	else if(bin_muedk>=1&&bin_nn==0) fit_params[bin_rec_type][bin_nn][bin_muedk][k]=fit_params[bin_rec_type][bin_nn][bin_muedk-1][k];
	else if(bin_nn>=1&&bin_muedk>=1) fit_params[bin_rec_type][bin_nn][bin_muedk][k]=(fit_params[bin_rec_type][bin_nn][bin_muedk-1][k]+fit_params[bin_rec_type][bin_nn-1][bin_muedk][k])/2.;
	else if(bin_muedk==0&&bin_nn>=1) fit_params[bin_rec_type][bin_nn][bin_muedk][k]=fit_params[bin_rec_type][bin_nn-1][bin_muedk][k];
	else                             fit_params[bin_rec_type][bin_nn][bin_muedk][k]=0;
      }
      }
    }
  }
  
  for(int bin_rec_type=0;bin_rec_type<nbins_rec_type;bin_rec_type++){    
    for(int bin_nn=0;bin_nn<nbins_nn;bin_nn++){    
      for(int bin_muedk=0;bin_muedk<nbins_muedk;bin_muedk++){
	delete bin_value[bin_rec_type][bin_nn][bin_muedk];
	delete g2[bin_rec_type][bin_nn][bin_muedk];
      }
    }
  }
  
  }

void Energy_recon::remake_tree(){
  Int_t nn;
  Int_t ipnu;
  Float_t amom;
  Float_t amom_forcorr;
  Float_t amom_pri;
  Float_t amom_corrected;
  Int_t itype;
  Int_t muedk;

  osc_tuple_rec=new TTree("osc_tuple_rec","osc_tuple_rec");
  osc_tuple_origin->SetBranchAddress("nn"     , &nn);
  osc_tuple_origin->SetBranchAddress("ipnu"   , &ipnu);
  osc_tuple_origin->SetBranchAddress("amom"   , &amom    );
  osc_tuple_origin->SetBranchAddress("amom_forcorr"   , &amom_forcorr    );
  osc_tuple_origin->SetBranchAddress("amom_pri"   , &amom_pri);
  osc_tuple_origin->SetBranchAddress("itype"  , &itype   );
  osc_tuple_origin->SetBranchAddress("muedk"  , &muedk   );


  
  osc_tuple_rec->Branch("amom_corrected"   , &amom_corrected    , "amom_corrected/F"    );
  
  Int_t nevent=osc_tuple_origin->GetEntries();
  Int_t table_type;
  Float_t mass;
  Int_t rec_type_bin;
  Int_t nnbin_input;
  Int_t muedkbin_input;
  Float_t fit_param_mean[5];
  Float_t nn_low;
  Float_t nn_high;
  Float_t muedk_low;
  Float_t muedk_high;
  for(Int_t i=0;i<nevent;i++){
    osc_tuple_origin->GetEntry(i);
    rec_type_bin=rec_type[(itype-1)%(EndOfTypes-1)];
    
    if(rec_type_bin==0||rec_type_bin==2) mass=0.511;
    if(rec_type_bin==1||rec_type_bin==3) mass=105.6;
    else mass=0;

    amom_corrected=0;
    if(rec_type_bin>=0&&rec_type_bin<=3){
      for(Int_t j=0;j<5;j++){
	if(nn<=6){
	  nnbin_input=nn;
	  nn_low=nn;
	  nn_high=nn+1;
	}      
	else{
	  nnbin_input=(Int_t)6.*(log((Float_t)nn/6.)/2.+1);
	  if(nnbin_input>=nbins_nn) nnbin_input=nbins_nn-1;
	  nn_low=6.*exp(2.*((Float_t)nnbin_input/6.-1));
	  nn_high=6.*exp(2.*((Float_t)(nnbin_input+1)/6.-1));
	}
	if(muedk<nbins_muedk) muedkbin_input=muedk;     
	else muedkbin_input=nbins_muedk-1;
	
	muedk_low=muedkbin_input;
	muedk_high=muedkbin_input+1;
	
	if(nnbin_input<nbins_nn-1&&muedkbin_input<nbins_muedk-1){
	  fit_param_mean[j]=(fit_params[rec_type_bin][nnbin_input][muedkbin_input][j]*(nn_high-nn)*(muedk_high-muedk)
			     +fit_params[rec_type_bin][nnbin_input+1][muedkbin_input][j]*(nn-nn_low)*(muedk_high-muedk)
			     +fit_params[rec_type_bin][nnbin_input][muedkbin_input+1][j]*(nn_high-nn)*(muedk-muedk_low)
			     +fit_params[rec_type_bin][nnbin_input+1][muedkbin_input+1][j]*(nn-nn_low)*(muedk-muedk_low))/(nn_high-nn_low)/(muedk_high-muedk_low);
	}
	else if(nnbin_input==nbins_nn-1&&muedkbin_input<nbins_muedk-1){
	  fit_param_mean[j]=(fit_params[rec_type_bin][nnbin_input][muedkbin_input][j]*(muedk_high-(muedk))
			     +fit_params[rec_type_bin][nnbin_input][muedkbin_input+1][j]*((muedk)-muedk_low))/(muedk_high-muedk_low);
	}
	else if(nnbin_input<nbins_nn-1&&muedkbin_input==nbins_muedk-1){
	  fit_param_mean[j]=(fit_params[rec_type_bin][nnbin_input][muedkbin_input][j]*(nn_high-nn)
			     +fit_params[rec_type_bin][nnbin_input+1][muedkbin_input][j]*(nn-nn_low))/(nn_high-nn_low);
	}
	else if(nnbin_input==nbins_nn-1&&muedkbin_input==nbins_muedk-1){
	  fit_param_mean[j]=fit_params[rec_type_bin][nnbin_input][muedkbin_input][j];
	}      
      }
      amom_corrected=(fit_param_mean[0]+fit_param_mean[1]*(amom_pri-pow(10,fit_param_mean[3]))
		      +(fit_param_mean[2]*(amom_pri>fit_param_mean[3])+fit_param_mean[4]*(amom_pri<fit_param_mean[3]))*pow((amom_pri-pow(10,fit_param_mean[3])),2));
      amom_corrected+=sqrt(amom_pri*amom_pri+mass*mass);
      if(rec_type_bin==2||rec_type_bin==3){
	amom_corrected+=(amom_forcorr-sqrt(amom_pri*amom_pri+mass*mass));
      }
      if(amom_corrected<amom_forcorr) amom_corrected=amom_forcorr;
    }
    else amom_corrected=0;
    osc_tuple_rec->Fill();
  }
}

void Energy_recon::offset_fix(Bool_t read_io){
  Int_t nn;
  Int_t ipnu;
  Float_t pnu;
  Int_t mode;
  Float_t amom_corrected;
  Float_t amom_corrected_tmp;
  Int_t itype;
  Int_t muedk;
  Double_t xq[1]={0.5};
  Double_t yq[1];
  Int_t type;
  Int_t target_type;
  Int_t rec_type_bin;
  Int_t bin_nn;
  Int_t bin_muedk;
  SetBranch_ERecon();

  osc_tuple_rec->SetBranchAddress("amom_corrected"   , &amom_corrected_tmp    );
  otree_erecon = new TTree("osc_tuple_erecon", "tree build for SK osc analyses");
  otree_erecon->Branch("amom_corrected"   , &amom_corrected, "amom_corrected");

  
  Int_t nevent=osc_tuple_origin->GetEntries();

  if(read_io==false){
  for(Int_t i=0;i<nbins_rec_type;i++){
    for(Int_t j=0;j<nbins_nn;j++){	
      for(Int_t k=0;k<nbins_muedk;k++){
	amom_offset_hist[i][j][k]=new TH1F(Form("amom_offset_hist_%d_%d_%d",k,j,i),"",200,-1,1);
      }
    }   
  }

  for(Int_t i=0;i<nevent;i++){
    osc_tuple_origin->GetEntry(i);
    osc_tuple_rec->GetEntry(i);
    nn=os2->nn;
    ipnu=os2->ipnu;
    pnu=os2->pnu;
    mode=os2->mode;
    itype=os2->itype;
    muedk=os2->muedk;
    rec_type_bin=rec_type[(itype-1)%(EndOfTypes-1)];    
    if(rec_type_bin==0||rec_type_bin==2)      target_type=12;
    else if(rec_type_bin==1||rec_type_bin==3) target_type=14;
    else continue;    

    if(nn<6) bin_nn=nn;
    else if(nn<nbins_nn) bin_nn=(Int_t)6*(log((Float_t)nn/6.)/2.+1);
    else continue;
    if(muedk<nbins_muedk) bin_muedk=muedk;
    else continue;

    if(pnu>Elow_resonance&&pnu<Ehigh_resonance&&abs(mode)<=27&&(abs(ipnu)==target_type)){
      amom_offset_hist[rec_type_bin][bin_nn][bin_muedk]->Fill((amom_corrected_tmp-pnu*1000.)/pnu/1000.);
    }    
  }
  
  for(int i=0;i<nbins_rec_type;i++){ 
    for(int j=0;j<nbins_nn;j++){
      for(int k=0;k<nbins_muedk;k++){
	if(amom_offset_hist[i][j][k]->Integral()>10){
	amom_offset_hist[i][j][k]->GetQuantiles(1,yq,xq);
	amom_offset[i][j][k]=yq[0];
	}
	else if(k>=1&&j==0) amom_offset[i][j][k]=amom_offset[i][j][k-1];	
	else if(j>=1&&k==0) amom_offset[i][j][k]=amom_offset[i][j-1][k];	
	else if(j>=1&&k>=1) amom_offset[i][j][k]=(amom_offset[i][j][k-1]+amom_offset[i][j-1][k])/2.;	
	else                amom_offset[i][j][k]=0;
      }
    }
  }
  }
  for(Int_t i=0;i<nevent;i++){
    osc_tuple_origin->GetEntry(i);
    osc_tuple_rec->GetEntry(i);
    nn=os2->nn;
    ipnu=os2->ipnu;
    pnu=os2->pnu;
    mode=os2->mode;
    itype=os2->itype;
    muedk=os2->muedk;
    rec_type_bin=rec_type[(os2->itype-1)%(EndOfTypes-1)];    

    if(nn<6) bin_nn=nn;
    else if(nn<nbins_nn) bin_nn=(Int_t)6*(log((Float_t)nn/6.)/2.+1);    
    else bin_nn=nbins_nn-1;
    
    if(muedk<nbins_muedk) bin_muedk=muedk;
    else bin_muedk=nbins_muedk-1;
    amom_corrected=amom_corrected_tmp/(1.+amom_offset[rec_type_bin][bin_nn][bin_muedk]);
    otree_erecon->Fill();
  }
}

Double_t Energy_recon::MAD(TH1F* h1, Double_t median_value){
  Int_t nbin=h1->GetXaxis()->GetNbins();
  Int_t nevent;
  if(h1->GetEntries()>0) nevent=h1->Integral(0,nbin);
  else nevent=0;
  Double_t minimum=h1->GetXaxis()->GetBinLowEdge(1);
  Double_t maximum=h1->GetXaxis()->GetBinLowEdge(nbin)+h1->GetXaxis()->GetBinWidth(nbin);
  Double_t value;
  Double_t xq[1]={0.5};
  Double_t yq[1];
  TH1F* h2 = new TH1F("h2","",nbin,minimum-median_value,maximum-median_value);
  for(Int_t i=0;i<nbin;i++){
    value=h1->GetBinContent(i+1);
    h2->Fill(fabs(h1->GetBinLowEdge(i+1)+h1->GetBinWidth(i+1)-median_value),value);
  }
  h2->GetQuantiles(1,yq,xq);
  h2->Delete();
  return yq[0];
}

void Energy_recon::SetBranch_ERecon(){

  osc_tuple_origin->SetBranchAddress("nn"     , &os2->nn);
  osc_tuple_origin->SetBranchAddress("nn_mctruth"     , &os2->nn_mctruth);
  osc_tuple_origin->SetBranchAddress("nring"  ,&os2->nring);

  osc_tuple_origin->SetBranchAddress("ipnu"   , &os2->ipnu);
  osc_tuple_origin->SetBranchAddress("dirnu"  , &os2->dirnu);
  osc_tuple_origin->SetBranchAddress("pnu"    , &os2->pnu);
  osc_tuple_origin->SetBranchAddress("mode"   , &os2->mode);
  osc_tuple_origin->SetBranchAddress("ip"     , &os2->ip);
  osc_tuple_origin->SetBranchAddress("dprob"  , &os2->dprob);
  osc_tuple_origin->SetBranchAddress("dir"    , &os2->dir);
  osc_tuple_origin->SetBranchAddress("amom"   , &os2->amom); 
  osc_tuple_origin->SetBranchAddress("amom_forcorr"   , &os2->amom_forcorr);
  osc_tuple_origin->SetBranchAddress("amom_pri"   , &os2->amom_pri);
  osc_tuple_origin->SetBranchAddress("path"   , &os2->path);
  osc_tuple_origin->SetBranchAddress("wall"   , &os2->wall);
  osc_tuple_origin->SetBranchAddress("itype"  , &os2->itype);
  osc_tuple_origin->SetBranchAddress("muedk"  , &os2->muedk);
  osc_tuple_origin->SetBranchAddress("flxg"   , &os2->flxg );
  osc_tuple_origin->SetBranchAddress("flxgo"  , &os2->flxgo);
  osc_tuple_origin->SetBranchAddress("flxh"   , &os2->flxh );
  osc_tuple_origin->SetBranchAddress("flxho"  , &os2->flxho);
  osc_tuple_origin->SetBranchAddress("weightx", &os2->weightx);
  osc_tuple_origin->SetBranchAddress("fsiweight", &os2->fsiweight);
  osc_tuple_origin->SetBranchAddress("oscweight3f", &os2->oscweight3f);
  osc_tuple_origin->SetBranchAddress("nn_selected", &os2->nn_selected);
  osc_tuple_origin->SetBranchAddress("nn_output", &os2->nn_output);  
  
}

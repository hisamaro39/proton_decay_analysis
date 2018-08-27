#include <math.h>

float calc_poisson(float rate,float exposure, float efficiency, float nbkg, int ncand);

void calc_limit(){

  float decay_rate;
  float exposure = 0.306;//Mton*years
  float nbkg = 0.87; 
  float nbkg2 = 5*0.87; 
  float eff = 0.2;
  int observed = 0;

  float poisson,poisson2,poisson_new;
  TGraph *g = new TGraph();
  TGraph *g2 = new TGraph();
  float allsum=0.,allsum2=0.;
  for(int i=0;i<10000;i++){
    decay_rate = i*0.01;
    float total_expected = exposure*decay_rate*eff+nbkg;
    float total_expected2 = exposure*decay_rate*eff+nbkg2;
    poisson = TMath::Poisson(observed,total_expected);
    poisson2 = TMath::Poisson(observed,total_expected2);
    poisson_new = calc_poisson(decay_rate,exposure,eff,nbkg,observed);
    cout << "poisson/new=" << poisson << "/" << poisson_new << endl;
    g->SetPoint(i,decay_rate,poisson);
    g2->SetPoint(i,decay_rate,poisson2);
    //cout << "decay_rate/total/poisson=" << decay_rate << "/" << total_expected << "/" << poisson << endl; 
    allsum += 0.1*poisson;
    allsum2 += 0.1*poisson2;
  }
  cout << "allsum/allsum2=" << allsum << "/" << allsum2 << endl;
  float sum=0.,sum2=0.;
  float rate_limit,rate_limit2;
  for(int i=0;i<10000;i++){
    decay_rate = i*0.01;//[(Mton*years)^-1]
    float total_expected = exposure*decay_rate*eff+nbkg;
    poisson = TMath::Poisson(observed,total_expected);
    sum += 0.1*poisson;
    if(sum/allsum > 0.9) {
      rate_limit=decay_rate;
      break;
    }
  }
  for(int i=0;i<10000;i++){
    decay_rate = i*0.01;
    float total_expected2 = exposure*decay_rate*eff+nbkg2;
    poisson2 = TMath::Poisson(observed,total_expected2);
    sum2 += 0.1*poisson2;
    if(sum2/allsum2 > 0.9) {
      rate_limit2=decay_rate;
      break;
    }
  }
  float mol_H2O  = exposure*1e12/18;
  float num_proton = mol_H2O * 6.02*1e23;
  cout << "mol_H2O=" << mol_H2O << endl;
  cout << "num_proton=" << num_proton << endl;
  cout << "rate_limit/rate_limit2=" << rate_limit << "/" << rate_limit2 << std::endl;
  float convert = num_proton/exposure;
  float life_time_limit = convert/rate_limit;
  float life_time_limit2 = convert/rate_limit2;
  cout << "life time limit is " << life_time_limit << " years" << endl;
  cout << "life time limit2 is " << life_time_limit2 << " years" << endl;

  //g->Draw("al");
  //g2->SetLineColor(2);
  //g2->Draw("l");


}

float calc_poisson(float rate,float exposure, float efficiency, float nbkg, int ncand){

    float m = (rate*exposure*efficiency+nbkg);

    float rcand = ncand * log(m) - math::lgamma(ncand + 1) ;

    float func = exp ( -m + rcand );

  return func;
}

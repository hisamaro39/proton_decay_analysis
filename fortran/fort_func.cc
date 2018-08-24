#include "fortran/fort_func.h"
#include <iostream>
#include <fstream>
#include <math.h>

#define PI 3.141592

fort_func::fort_func(){
  std::cout << "fort_func::Initialize" << std::endl;

  std::ifstream table("/disk01/usr5/matanaka/proton_decay/Analysis/fortran/hkkm11mt.dat");
  std::ifstream tablelow("/disk01/usr5/matanaka/proton_decay/Analysis/fortran/hkkm11low.dat");

  std::string line,linelow;
  for(int s=0;s<2;s++){//solar activity
    for(int t=0;t<20;t++){//zenith angle
      for(int p=0;p<12;p++){//azmithal angle
        std::getline(table, line);//just for unneed line
        std::getline(table, line);//just for unneed line
        std::getline(tablelow, linelow);//just for unneed line
        std::getline(tablelow, linelow);//just for unneed line
        for(int e=0;e<101;e++){//energy
          std::getline(table, line);
          std::getline(tablelow, linelow);
          sscanf(line.c_str(), "%lf %lf %lf %lf %lf", &energy_range[e], &bfnu[s][t][p][e][0], &bfnu[s][t][p][e][1], &bfnu[s][t][p][e][2], &bfnu[s][t][p][e][3]);
          sscanf(linelow.c_str(), "%lf %lf %lf %lf %lf", &energy_range_low[e], &bfnu_low[s][t][p][e][0], &bfnu_low[s][t][p][e][1], &bfnu_low[s][t][p][e][2], &bfnu_low[s][t][p][e][3]);
        }
      }
    }
  }

  std::ifstream table2("/disk01/usr5/matanaka/proton_decay/Analysis/fortran/honda3d_prodhgt");
  std::string line2;
  double energy,value;
  for(int f=0;f<4;f++){//flavor
    for(int z=0;z<20;z++){//zenith angle
      for(int e=0;e<31;e++){//energy
        table2 >> energy;
        //std::cout << "energy=" << energy << std::endl;
        energy_range2[e] = energy;
        for(int v=0;v<19;v++){
          table2 >> value;
          //std::cout << "v/e/z/f=" << v << "/" << e << "/" << z << "/" << f << std::endl;
          //std::cout << "value=" << value << std::endl;
          h3d[v][e][z][f] = value;
        }
      }
    }
  }

}


double fort_func::calc_flux(float d_pnu,float *d_dirnu,float Solar,int nu){
  //std::cout << "fort_func::calc_flux !!" << std::endl;

  //std::cout << "d_pnu/Sloar/nu=" << d_pnu << "/" << Solar << "/" << nu << std::endl;
  //std::cout << "d_dirnu 0/1/2=" << d_dirnu[0] << "/" << d_dirnu[1] << "/" << d_dirnu[2] << std::endl;

  int ipkind = -1;
  if(nu==12) ipkind = 0;
  if(nu==-12) ipkind = 1;
  if(nu==14) ipkind = 2;
  if(nu==-14) ipkind = 3;

  float new_dirnu[3];
  new_dirnu[2] = -1*d_dirnu[2];
  new_dirnu[0] = -1*(cos(PI*2435./60./180.)*d_dirnu[0] +sin(PI*2435./60./180.)*d_dirnu[1]);
  new_dirnu[1] = sin(PI*2435./60./180.)*d_dirnu[0] - cos(PI*2435./60./180.)*d_dirnu[1];

  if(d_pnu > energy_range[100]) return calc_flux_high(d_pnu, new_dirnu, Solar, ipkind);
  if(d_pnu < energy_range[0]) return calc_flux_low(d_pnu, new_dirnu, Solar, ipkind);
  else return calc_flux_middle(d_pnu, new_dirnu, Solar, ipkind);


}

double fort_func::calc_flux_low(float d_pnu,float *new_dirnu,float Solar,int ipkind){
  //std::cout << "fort_func::calc_flux_low !!" << std::endl;

  int ebin;
  for(ebin=0;ebin<101;ebin++)
    if(d_pnu<energy_range_low[ebin]) break;
  ebin = ebin-1;

  int zenith_bin_high = (int) (new_dirnu[2]*10 + 10.5);
  int zenith_bin_low = zenith_bin_high - 1;
  //std::cout << "new_dir[2]=" << new_dirnu[2] << std::endl;
  if(zenith_bin_low==-1) zenith_bin_low = 0;
  if(new_dirnu[2]<-0.95){
    zenith_bin_low = 0;
    zenith_bin_high = 0;
  }
  if(new_dirnu[2]>0.95){
    zenith_bin_low = 19;
    zenith_bin_high = 19;
  }

  float phi = -99;
  if(fabs(new_dirnu[0])<0.0 && fabs(new_dirnu[1])<0.0) phi = 0;
  else phi = atan2(new_dirnu[1],new_dirnu[0]);
  //std::cout << "phi=" << phi << std::endl;

  float new_phi=-99;
  if(phi<0) new_phi = phi + 2*PI; 
  else new_phi = phi; 
  //std::cout << "new_phi=" << new_phi << std::endl;

  int phibin_low= (int) (new_phi/(PI/6)-0.5) ;
  int phibin_high = phibin_low + 1;
  if(fabs(phi)<PI/12){
    phibin_low = 11;
    phibin_high = 0;
  }

  float zenith_low = (zenith_bin_low+1)*0.1 - 1.05;
  float zenith_high = (zenith_bin_high+1)*0.1 - 1.05;

  //std::cout << "zenithbin low/high=" << zenith_bin_low << "/" << zenith_bin_high << std::endl;
  //std::cout << "zenith low/high=" << zenith_low << "/" << zenith_high << std::endl;

  float phi_low = new_phi - ( (PI/6)*(phibin_low - 0.5) );
  float phi_high = (PI/6)*(phibin_high - 0.5) - new_phi;
  if(phibin_low==11){
    phi_low = phi + PI/12;
    phi_high = PI/12 - phi;
  }
  //std::cout << "phibin low/high=" << phibin_low << "/" << phibin_high << std::endl;
  //std::cout << "phi low/high=" << phi_low << "/" << phi_high << std::endl;

  float fxleltpl = bfnu_low[0][zenith_bin_low][phibin_low][ebin][ipkind]*(1-Solar) 
    + bfnu_low[1][zenith_bin_low][phibin_low][ebin][ipkind]*Solar; 
  float fxleltph = bfnu_low[0][zenith_bin_low][phibin_high][ebin][ipkind]*(1-Solar) 
    + bfnu_low[1][zenith_bin_low][phibin_high][ebin][ipkind]*Solar; 
  float flxlelt = (fxleltpl-fxleltph)/(PI/6)*phi_low + fxleltpl;
  //std::cout << "fxleltpl=" << fxleltpl << std::endl;
  //std::cout << "fxleltph=" << fxleltph << std::endl;
  //std::cout << "flxlelt=" << flxlelt << std::endl;
  float fxlehtpl = bfnu_low[0][zenith_bin_high][phibin_low][ebin][ipkind]*(1-Solar) 
    + bfnu_low[1][zenith_bin_high][phibin_low][ebin][ipkind]*Solar; 
  float fxlehtph = bfnu_low[0][zenith_bin_high][phibin_high][ebin][ipkind]*(1-Solar) 
    + bfnu_low[1][zenith_bin_high][phibin_high][ebin][ipkind]*Solar; 
  float flxleht = (fxlehtpl-fxlehtph)/(PI/6)*phi_low + fxlehtpl;

  float flxle = (flxleht-flxlelt)/(zenith_high-zenith_low)*(new_dirnu[2]-zenith_low)+flxlelt;

  float fxheltpl = bfnu_low[0][zenith_bin_low][phibin_low][ebin+1][ipkind]*(1-Solar) 
    + bfnu_low[1][zenith_bin_low][phibin_low][ebin+1][ipkind]*Solar; 
  float fxheltph = bfnu_low[0][zenith_bin_low][phibin_high][ebin+1][ipkind]*(1-Solar) 
    + bfnu_low[1][zenith_bin_low][phibin_high][ebin+1][ipkind]*Solar; 
  float flxhelt = (fxheltpl-fxheltph)/(PI/6)*phi_low + fxheltpl;

  float fxhehtpl = bfnu_low[0][zenith_bin_high][phibin_low][ebin+1][ipkind]*(1-Solar) 
    + bfnu_low[1][zenith_bin_high][phibin_low][ebin+1][ipkind]*Solar; 
  float fxhehtph = bfnu_low[0][zenith_bin_high][phibin_high][ebin+1][ipkind]*(1-Solar) 
    + bfnu_low[1][zenith_bin_high][phibin_high][ebin+1][ipkind]*Solar; 
  float flxheht = (fxhehtpl-fxhehtph)/(PI/6)*phi_low + fxhehtpl;

  float flxhe = (flxheht-flxhelt)/(zenith_high-zenith_low)*(new_dirnu[2]-zenith_low)+flxhelt;

  if(zenith_bin_low == zenith_bin_high){
    flxle = flxlelt;
    flxhe = flxhelt;
  }

  float B = log(flxle/flxhe) / log(energy_range_low[ebin]/energy_range_low[ebin+1]);
  float A = flxle/pow(energy_range_low[ebin],B);

  float ans = A*pow(d_pnu,B); 

  return ans;
}

double fort_func::calc_flux_middle(float d_pnu,float *new_dirnu,float Solar,int ipkind){
  //std::cout << "fort_func::calc_flux_middle !!" << std::endl;
  int ebin;
  for(ebin=0;ebin<101;ebin++)
    if(d_pnu<energy_range[ebin]) break;
  ebin = ebin-1;

  int zenith_bin_high = (int) (new_dirnu[2]*10 + 10.5);
  int zenith_bin_low = zenith_bin_high - 1;
  //std::cout << "new_dir[2]=" << new_dirnu[2] << std::endl;
  if(zenith_bin_low==-1) zenith_bin_low = 0;
  if(new_dirnu[2]<-0.95){
    zenith_bin_low = 0;
    zenith_bin_high = 0;
  }
  if(new_dirnu[2]>0.95){
    zenith_bin_low = 19;
    zenith_bin_high = 19;
  }

  float phi = -99;
  if(fabs(new_dirnu[0])<0.0 && fabs(new_dirnu[1])<0.0) phi = 0;
  else phi = atan2(new_dirnu[1],new_dirnu[0]);
  //std::cout << "phi=" << phi << std::endl;

  float new_phi=-99;
  if(phi<0) new_phi = phi + 2*PI; 
  else new_phi = phi; 
  //std::cout << "new_phi=" << new_phi << std::endl;

  int phibin_low= (int) (new_phi/(PI/6)-0.5) ;
  int phibin_high = phibin_low + 1;
  if(fabs(phi)<PI/12){
    phibin_low = 11;
    phibin_high = 0;
  }

  float zenith_low = (zenith_bin_low+1)*0.1 - 1.05;
  float zenith_high = (zenith_bin_high+1)*0.1 - 1.05;

  //std::cout << "zenithbin low/high=" << zenith_bin_low << "/" << zenith_bin_high << std::endl;
  //std::cout << "zenith low/high=" << zenith_low << "/" << zenith_high << std::endl;

  float phi_low = new_phi - ( (PI/6)*(phibin_low - 0.5) );
  float phi_high = (PI/6)*(phibin_high - 0.5) - new_phi;
  if(phibin_low==11){
    phi_low = phi + PI/12;
    phi_high = PI/12 - phi;
  }
  //std::cout << "phibin low/high=" << phibin_low << "/" << phibin_high << std::endl;
  //std::cout << "phi low/high=" << phi_low << "/" << phi_high << std::endl;

  float fxleltpl = bfnu[0][zenith_bin_low][phibin_low][ebin][ipkind]*(1-Solar) 
    + bfnu[1][zenith_bin_low][phibin_low][ebin][ipkind]*Solar; 
  float fxleltph = bfnu[0][zenith_bin_low][phibin_high][ebin][ipkind]*(1-Solar) 
    + bfnu[1][zenith_bin_low][phibin_high][ebin][ipkind]*Solar; 
  float flxlelt = (fxleltpl-fxleltph)/(PI/6)*phi_low + fxleltpl;
  //std::cout << "fxleltpl=" << fxleltpl << std::endl;
  //std::cout << "fxleltph=" << fxleltph << std::endl;
  //std::cout << "flxlelt=" << flxlelt << std::endl;
  float fxlehtpl = bfnu[0][zenith_bin_high][phibin_low][ebin][ipkind]*(1-Solar) 
    + bfnu[1][zenith_bin_high][phibin_low][ebin][ipkind]*Solar; 
  float fxlehtph = bfnu[0][zenith_bin_high][phibin_high][ebin][ipkind]*(1-Solar) 
    + bfnu[1][zenith_bin_high][phibin_high][ebin][ipkind]*Solar; 
  float flxleht = (fxlehtpl-fxlehtph)/(PI/6)*phi_low + fxlehtpl;

  float flxle = (flxleht-flxlelt)/(zenith_high-zenith_low)*(new_dirnu[2]-zenith_low)+flxlelt;

  float fxheltpl = bfnu[0][zenith_bin_low][phibin_low][ebin+1][ipkind]*(1-Solar) 
    + bfnu[1][zenith_bin_low][phibin_low][ebin+1][ipkind]*Solar; 
  float fxheltph = bfnu[0][zenith_bin_low][phibin_high][ebin+1][ipkind]*(1-Solar) 
    + bfnu[1][zenith_bin_low][phibin_high][ebin+1][ipkind]*Solar; 
  float flxhelt = (fxheltpl-fxheltph)/(PI/6)*phi_low + fxheltpl;

  float fxhehtpl = bfnu[0][zenith_bin_high][phibin_low][ebin+1][ipkind]*(1-Solar) 
    + bfnu[1][zenith_bin_high][phibin_low][ebin+1][ipkind]*Solar; 
  float fxhehtph = bfnu[0][zenith_bin_high][phibin_high][ebin+1][ipkind]*(1-Solar) 
    + bfnu[1][zenith_bin_high][phibin_high][ebin+1][ipkind]*Solar; 
  float flxheht = (fxhehtpl-fxhehtph)/(PI/6)*phi_low + fxhehtpl;

  float flxhe = (flxheht-flxhelt)/(zenith_high-zenith_low)*(new_dirnu[2]-zenith_low)+flxhelt;

  if(zenith_bin_low == zenith_bin_high){
    flxle = flxlelt;
    flxhe = flxhelt;
  }

  float B = log(flxle/flxhe) / log(energy_range[ebin]/energy_range[ebin+1]);
  float A = flxle/pow(energy_range[ebin],B);

  float ans = A*pow(d_pnu,B); 

  return ans;
}

double fort_func::calc_flux_high(float d_pnu,float *new_dirnu,float Solar,int ipkind){
  //std::cout << "fort_func::calc_flux_high !!" << std::endl;

  return 0.;
}

void  fort_func::nebaseline(float *path, int nutype, float costheta, float E_neu, float d){
  int ii;
  float R = 6371.000;//Radius of the Earth in km
  float a = 1;
  float x = R+d;
  float b = (costheta < 0)? 2.*x*costheta : -2.*x*costheta;
  float c = x*x - R*R;
  float Q=0.;//Distance to surface along zenith angle
  if(b*b < 4*a*c) Q = sqrt(2*d*x);//only true very near horizon, use tangent value there
  else if (costheta > 0) Q = (-1.*b - sqrt(b*b-4*a*c))/(2*a);
  else  Q = (-1.*b + sqrt(b*b-4*a*c))/(2*a);
  float cosalpha = fabs( (Q*Q+R*R-x*x)/(2.*Q*R) );//Angle at surface
  float hgt[20];
  neprodhgt_h3d(hgt,nutype,costheta,E_neu);

  float P=0.;
  if(fabs(costheta) < sqrt( (x*x-R*R)/(x*x) ) ) {
    for(int ii=0;ii<20;ii++)
      path[ii] = -1.*x*costheta + sqrt(pow(x*costheta,2)+pow(R+hgt[ii],2)-x*x);
  }
  else{
    for(int ii=0;ii<20;ii++){
      P = sqrt(pow(R*cosalpha,2)+pow(R+hgt[ii],2)-R*R)-R*cosalpha;
      if (costheta < 0.) 
        path[ii] = Q+P;
      else 
        path[ii] = P-Q;
    }
  }


}

void fort_func::neprodhgt_h3d(float *hgt,int nutype,float costheta,float E_nu){
  int costag = (int) ( (costheta+1)*10 ) + 1;
  if(costag > 20) costag = 20;
  float delta_cos = costheta - ( (costag-11.)/10 + 0.05 );

  int enetag = (int) (10*log10(E_nu) + 11);
  if(enetag < 0) enetag = 1;
  if(enetag > 31) enetag = 31;
  float de = pow(10,((enetag+1)-11.0)/10.0) - pow(10,(enetag-11.0)/10.0);
  float delta_e =  E_nu - pow(10,(enetag-11.0)/10.0);

  int nutag=-1;
  if(nutype==1) nutag = 2;
  if(nutype==2) nutag = 0;
  if(nutype==3 || nutype==-1) nutag = 3;
  if(nutype==4 || nutype==-2) nutag = 1;

  for(int ii=0;ii<19;ii++){
    float delta_h=0.,dhe=0.,dhc=0.;

    if(enetag < 31 && delta_e > 0.0) {
      dhe = h3d[ii][enetag][costag-1][nutag] - h3d[ii][enetag-1][costag-1][nutag];
      delta_h = delta_h + delta_e * dhe/de;
    }

    if( (costag==1 && delta_cos<0.) || (costag==20 && delta_cos>0.) ) dhc = 0.;
    else {
      if(delta_cos>0.) dhc = h3d[ii][enetag-1][costag][nutag] - h3d[ii][enetag-1][costag-1][nutag];
      else dhc = h3d[ii][enetag-1][costag-2][nutag] - h3d[ii][enetag-1][costag-1][nutag];
      delta_h = delta_h + fabs(delta_cos) * dhc/0.1;
    }

    hgt[ii] = h3d[ii][enetag-1][costag-1][nutag] + delta_h;

  }
  hgt[19] = hgt[18];
}





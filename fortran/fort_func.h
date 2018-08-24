#ifndef __fort_func__
#define __fort_func__

#include <iostream>
#include <string>

class fort_func
{
 public :

   fort_func();

   double calc_flux(float d_pnu,float *d_dirnu,float Solar,int nu);
   double calc_flux_low(float d_pnu,float *new_dirnu,float Solar,int ipkind);
   double calc_flux_middle(float d_pnu,float *new_dirnu,float Solar,int ipkind);
   double calc_flux_high(float d_pnu,float *new_dirnu,float Solar,int ipkind);
   void nebaseline(float *path, int nutype, float costheta, float E_neu, float d);
   void neprodhgt_h3d(float *hgt,int nutype,float costheta,float E_neu);

   double energy_range[101];
   double energy_range_low[101];
   double energy_range2[31];
   double bfnu[2][20][12][101][4]; 
   double bfnu_low[2][20][12][101][4]; 
   double h3d[19][31][20][4]; 

 protected: 

};

#endif

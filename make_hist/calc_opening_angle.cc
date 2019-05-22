void calc_opening_angle(){

  gROOT->SetStyle("Plain");
  gStyle->SetOptStat(0);
  gStyle->SetPadBottomMargin(0.1);

  float mom_e = 39.7;
  float mom_m = 213.7;

  float n = 1.333;
  float mass_e = 0.511;//MeV
  float mass_m = 105.7;//MeV
  float theta_m = acos(sqrt(1+pow(mass_m,2)/pow(mom_m,2))/n);
  float angle_m = 180*theta_m/3.14159;
  float theta_e = acos(sqrt(1+pow(mass_e,2)/pow(mom_e,2))/n);
  float angle_e = 180*theta_e/3.14159;
  
  cout << "electron mom/angle=" << "/" << mom_e << "/" << angle_e << endl;
  cout << "muon mom/angle=" << "/" << mom_m << "/" << angle_m << endl;
 
}

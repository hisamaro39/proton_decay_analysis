void syst_summary(){
  syst_signal[6][4][2][9];//mode,period,box,type
  syst_bkg[4][4][10];//mode,period,type

  //Signal
  //p->eee sk1 low
  syst_signal[0][0][0][0] = 4.0;//Correlated Decay
  syst_signal[0][0][0][1] = 10.1;//Fermi Motion
  syst_signal[0][0][0][2] = 0.3;//Fiductial Volume
  syst_signal[0][0][0][3] = 1.5;//Detector Non-Uniformity
  syst_signal[0][0][0][4] = 4.5;//Energy Scale
  syst_signal[0][0][0][5] = 0.5;//PID
  syst_signal[0][0][0][6] = 2.5;//Ring Counting
  syst_signal[0][0][0][7] = 2.0;//Decay Electron
  syst_signal[0][0][0][8] = 0;//Neutron Tagging


}

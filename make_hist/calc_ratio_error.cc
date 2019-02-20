void calc_ratio_error(){

  float value1[] = {0.45,0.22,0.11,0.14};
  float value1_err[] = {0.06,0.03,.02,0.05};
  float value2[] = {0.09,0.04,0.03,0.06};
  float value2_err[] = {0.03,0.01,0.01,0.03};
  int nLoop = sizeof(value1)/sizeof(float);
  for (int n=0;n<nLoop;n++){
    float diff = value2[n]/value1[n] - 1;
    float diff_err = sqrt(pow(value2[n]*value1_err[n],2)+pow(value1[n]*value2_err[n],2))/pow(value1[n],2);
    cout << "value1=" << value1[n] << "+-" << value1_err[n] << " value2=" << value2[n] << "+-" << value2_err[n] << endl; 
    cout << "diff=" << diff << "+-" << diff_err << endl;
  }

}

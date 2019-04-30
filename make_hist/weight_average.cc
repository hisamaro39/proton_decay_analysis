void weight_average(){

  float value[3] = {13.5,2.9,22.3};
  float error[3] =  {2.6,0.7,10.0};
  float bunsi=0,bunbo=0;
  for(int i=0;i<3;i++){
    bunsi += value[i]*pow(1./error[i],2);
    bunbo += pow(1./error[i],2);
  }
  float ans = bunsi/bunbo;
  cout << "average is " << ans << endl;

}

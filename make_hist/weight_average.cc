void weight_average(){

  float value[4] = {3.3,2.9,0.9,2.6};
  float error[4] =  {0.8,0.8,0.3,0.7};
  float bunsi=0,bunbo=0;
  for(int i=0;i<4;i++){
    bunsi += value[i]*pow(1./error[i],2);
    bunbo += pow(1./error[i],2);
  }
  float ans = bunsi/bunbo;
  cout << "average is " << ans << endl;

}

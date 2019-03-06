void test_poisson(){
  float exp_bkg = 0.03;
  cout << "expected background is " << exp_bkg << endl;
  
  TGraph *g = new TGraph();
  TGraph *g2 = new TGraph();
  TH1 *h_poisson = new TH1F("h_poisson",";;",100,0,10);
  float sum=0;
  for(int i=0;i<10000;i++){
    float value = 0.01*i;
    float pois = TMath::Poisson(value,exp_bkg);
    g->SetPoint(i,value,pois);
    h_poisson->Fill(value,pois);
    sum += 0.01*pois;
  }
  //g->Draw("al");
  //h_poisson->Draw();
  //cout << "0 value = " << TMath::Poisson(exp_bkg,0) << "/" << TMath::PoissonI(exp_bkg,0) << endl;
  //cout << "sum = " << sum << endl;

  //calc integration
  float sum2=0;
  for(int j=0;j<10000;j++){
    float value = 0.01*j;
    float pois = TMath::Poisson(value,exp_bkg);
    g->SetPoint(j,value,pois);
    h_poisson->Fill(value,pois);
    sum2 += 0.01*pois;
    float ratio = sum2/sum;
    if(ratio>0.9999994){//for 5sigma discovery
      cout << "ratio/value=" << ratio << "/" << value << endl;
      break;
    }
  }
}

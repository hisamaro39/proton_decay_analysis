void test_smear(){
  gROOT->SetStyle("Plain");
  gStyle->SetOptStat(0);

  TRandom *generator = new TRandom();
  generator->SetSeed(0);
  float reso_vertex_r = gRandom->Gaus(0,1);
  TH1 *func = new TH1F("func",";;",120,-10,110);
  TH1 *func_smear = new TH1F("func_smear",";;",120,-10,110);
  TF1 *f1 = new TF1("f1","x*x",0,100);
  for(int i=0;i<1000000;i++){
    //float x = 100*generator->Rndm();
    //float x = generator->Exp(2);
    float x = f1->GetRandom();
    //float x_reso = gRandom->Gaus(0,10);
    float x_reso = (gRandom->Rndm()<0.5)? -1 : 1;
    //while(x_smear<0){
      x_smear = x + x_reso;
    //}
    func->Fill(x);
    func_smear->Fill(x_smear);
  }
  //TH1* frame = new TH1F();
  //frame->DrawFrame(-10, -1000, 110, 11000);
  func->Draw();
  func_smear->SetLineColor(2);
  func_smear->Draw("same");
  TH1 *ratio = (TH1*) func->Clone("func_clone");
  ratio->Divide(func_smear);
  //ratio->Draw();
}

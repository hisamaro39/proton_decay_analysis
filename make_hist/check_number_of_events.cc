void check_number_of_events(){

  TChain *fcmc = new TChain("h1");
  TChain *fcmc_apr16 = new TChain("h1");
  TChain *fcmc_may19 = new TChain("h1");
  TChain *syst_friend = new TChain("oscfitted_weights");
  TChain *fsi_friend = new TChain("fsifriend");

  for(int a=0;a<500;a++){
    if(a<10){
      fcmc->Add(Form("/disk02/atmpd6/sk4_dst/may19/fc_apr16mc/ntuple/000%d.sk4.presel.apfit.mccomb.ntag.root",a));
      fcmc_apr16->Add(Form("/disk01/atmpd5/sk4_dst/apr16/fc_mc/ntuple/mar16sk4.reduc.000%d_fQv5r0.root",a));
      //syst_friend->Add(Form("/disk01/usr3/raw/shared/fitted.error.friends/16b.201606.skmeeting/sk4/fcmc/mar16sk4.reduc.000%d_fQv5r0_ntag16c.osc.fitted.friend.root",a));
      syst_friend->Add(Form("/disk01/usr3/raw/shared/fitted.error.friends/repro.16b.fitted.error.friends/sk4/fcmc/000%d.osc.fitted.friend.root.presel.apfit.mccomb.ntag.root",a));
    }
    else if(a<100) {
      fcmc->Add(Form("/disk02/atmpd6/sk4_dst/may19/fc_apr16mc/ntuple/00%d.sk4.presel.apfit.mccomb.ntag.root",a));
      fcmc_apr16->Add(Form("/disk01/atmpd5/sk4_dst/apr16/fc_mc/ntuple/mar16sk4.reduc.00%d_fQv5r0.root",a));
      //syst_friend->Add(Form("/disk01/usr3/raw/shared/fitted.error.friends/16b.201606.skmeeting/sk4/fcmc/mar16sk4.reduc.00%d_fQv5r0_ntag16c.osc.fitted.friend.root",a));
      syst_friend->Add(Form("/disk01/usr3/raw/shared/fitted.error.friends/repro.16b.fitted.error.friends/sk4/fcmc/00%d.osc.fitted.friend.root.presel.apfit.mccomb.ntag.root",a));
    }
    else if(a<1000) {
      fcmc->Add(Form("/disk02/atmpd6/sk4_dst/may19/fc_apr16mc/ntuple/0%d.sk4.presel.apfit.mccomb.ntag.root",a));
      fcmc_apr16->Add(Form("/disk01/atmpd5/sk4_dst/apr16/fc_mc/ntuple/mar16sk4.reduc.0%d_fQv5r0.root",a));
      //syst_friend->Add(Form("/disk01/usr3/raw/shared/fitted.error.friends/16b.201606.skmeeting/sk4/fcmc/mar16sk4.reduc.0%d_fQv5r0_ntag16c.osc.fitted.friend.root",a));
      syst_friend->Add(Form("/disk01/usr3/raw/shared/fitted.error.friends/repro.16b.fitted.error.friends/sk4/fcmc/0%d.osc.fitted.friend.root.presel.apfit.mccomb.ntag.root",a));
    }
  }
  for(int b=2500;b<3000;b++){
    //fsi_friend->Add(Form("/disk01/usr5/raw/fsi.friends/16b_fsi_sets/sk4/fcmc/apr18_sk4_apfit_fcred_ntag_%d.fsifriend.root",b));
    fsi_friend->Add(Form("/home/matanaka/disk02_usr6/fsi/fsi_friend/sk4/file/fsi_friend.%d.root",b));
    fcmc_may19->Add(Form("/disk02/atmpd6/sk4_dst/may19/fc_mc/ntuple/%d.photon.fcred.ntag.root",b));
  }
  cout << "# of events in fcmc is " << fcmc->GetEntries() << endl;
  cout << "# of events in fcmc_apr16 is " << fcmc_apr16->GetEntries() << endl;
  cout << "# of events in fcmc_may19 is " << fcmc_may19->GetEntries() << endl;
  cout << "# of events in syst_friend is " << syst_friend->GetEntries() << endl;
  cout << "# of events in fsi_friend is " << fsi_friend->GetEntries() << endl;


}

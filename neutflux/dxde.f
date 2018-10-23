      Real*8 Function Dxde(Emu)
C      Subroutine Dxde(Emu,RDedx)
C                
C     Calculate Muon Dx/De IN The Standard Rock
C     [Data: W.Lohman, R.Kopp and R.Voss, Preprint CERN-85-03]         
C                
C     Input : EMU   ; Muon energy [GeV]     
C     Output: RANGE ; Muon range [g/cm^2]            
C                
      Implicit Real*8 (A-H,O-Z)             
      Dimension Emus(61),Dedx(61),Ccc(61)
      Data Emus/0.116,0.126,0.136,0.146,0.156,0.166,0.176,0.186,       
     &   0.2,  0.3,  0.4,  0.5,  0.6,  0.7,  0.8,  0.9,       
     &    1.,   2.,   4.,   6.,   8.,  10.,  15.,  20.,       
     &   30.,  40.,  50.,  60.,  70.,  80.,  90., 100.,       
     &  120., 140., 160., 180., 200., 220., 240., 260.,       
     &  280., 300., 350., 400., 450., 500., 600., 700.,       
     &  800., 900.,1000.,1500.,2000.,3000.,4000.,5000.,       
     & 6000.,7000.,8000.,9000.,10000./      
      Data Dedx/6.457,4.013,3.137,2.690,2.422,2.244,2.118,2.026,       
     & 1.933,1.713,1.687,1.696,1.714,1.733,1.752,1.771,       
     & 1.788,1.907,2.024,2.090,2.135,2.170,2.234,2.282,       
     & 2.355,2.416,2.470,2.520,2.569,2.615,2.661,2.706,       
     & 2.794,2.881,2.966,3.052,3.137,3.222,3.306,3.391,       
     & 3.475,3.560,3.771,3.983,4.194,4.406,4.831,5.258,       
     & 5.685,6.113,6.543,8.702,10.88,15.25,19.66,24.08,       
     & 28.52,32.97,37.43,41.90,46.37/       
      Data Aaa/0.00167D0/            
      Data Bbb/4.47D-6/            
      Data Elimit/10000.D0/           
      Data Nnn/61/
      Data Yp1/1.D30/
      Data Ypn/1.D30/       
      Data Ifl/0/         
C                
C----- Initialization ----------------      
      If(Ifl.EQ.0) Then
        Ifl=1
C --- Numerical Recipes
        Call Dspline(Emus,Dedx,Nnn,Yp1,Ypn,Ccc)
      End If

      If(Emu.LT.0.116D0) Then
        Dxde=0.D0
        Return
      End If
C                
c      If(Emu.GE.Elimit) Then
c        RDedx=(Aaa+Bbb*Emu)*1.D3
c      Else
        Call Dsplint(Emus,Dedx,Ccc,Nnn,Emu,RDedx)
c      End If
      Dxde=1.D0/RDedx

      Return 
      End     

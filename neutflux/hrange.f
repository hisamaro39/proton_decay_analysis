C   24/09/91 209141954  MEMBER NAME  HRANGE   *.FORT     M  E2FORT     
      Subroutine Hrange(Emu,Eth,Range)      
C                
C     CALCULATE MUON RANGE IN THE STANDARD ROCK      
C     [Data: W.Lohman, R.Kopp and R.Voss, Preprint CERN-85-03]         
C                
C     Originated by Y. Oyama       
C     Modified by M. Mori   24-SEP-1991 [SSL2 => NUMPAC]      
C        M. Mori   30-JAN-1992 [return 0 if EMU<ETH]          
C                
C     Input : EMU   ; Muon energy [GeV]     
C     ETH   ; Muon energy threshold [GeV]    
C     Output: RANGE ; Muon range [g/cm^2]            
C                
      Implicit Real*8 (A-H,O-Z)             
      External Dxde

      Data Emin /0.116D0/
      Data Emax /10000.D0/
      DATA Aaa/0.00167D0/            
      DATA Bbb/4.47E-6/
C                
      Eee=Aaa/Bbb

      If(Eth.LT.0.116) Then
        Range = 0.D0      
        Return
      End If 
C                
      If(Emu.LT.Eth) Then
        Range = 0.D0      
        Return
      End If               
C                
      If(Emu.GE.Emax) Then
        Call Dqgaus(Dxde,Emin,Emax,Smu)
        Smu=Smu*1.D3+(Dlog(1.+Emu/Eee)-Dlog(1.+Emax/Eee))/Bbb
      Else
        Call Dqgaus(Dxde,Emin,Emu,Smu)
        Smu=Smu*1.D3
      End If

      If(Eth.GE.Emax) Then
        Call Dqgaus(Dxde,Emin,Emax,Sth)
        Sth=Sth*1.D3+(Dlog(1.+Eth/Eee)-Dlog(1.+Emax/Eee))/Bbb
      Else
        Call Dqgaus(Dxde,Emin,Eth,Sth)
        Sth=Sth*1.D3
      End If

C      Range=(Smu-Sth)*1.D3
      Range=Smu-Sth
C      Range=Smu
C      write(1,*) smu,sth
      Return 
      End     







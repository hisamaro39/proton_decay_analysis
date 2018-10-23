************************************************************************
*     ---------------------------------
      FUNCTION RNNEUE(IPAR)
*     ---------------------------------
C
*     (Purpose)
C       GIVE RANDOM ENERGY OF NEUTRINO
C       DISTRIBUTION IS DUE TO FLUX*CROSS SECTION
C
*     (Input)
*       IPAR  ; KIND OF PARTICLE
*
C     (Output)
C       RNNEUE
C
*     (Creation Date and Author)
C       1983-11-22 ; M.NAKAHATA
C       1985-01-12 ; M.NAKAHATA
C                     PAULI PRINCIPLE AND FERMI MOTION IS CONSIDERED.
C                     USEABLE FOR  10MEV < NEU E  < 5GEV
C                     FOR 120MEV < NEU MU < 5GEV
C       1985-11-09 ; M.NAKAHATA
C                     ENEU > 5GEV WAS SUPPORTED
C       1988-06-25 ; T.KAJITA 
C                     RLOW, RMID WAS CHANGED IN ORDER TO CONSIDER
C                     SOLAR ACTIVITY AUTOMATICALLY.
C                     SOLAR MIN,MAX OR ANY COMBINATION OF MIN/MAX IS
C                     AVAILABLE
C                       IF(SOLACT=1.) -- SOLAR MAX.
C                       IF(SOLACT=0.) -- SOLAR MIN.
C                     ( IF(SOLACT=0.5) -- AVERAGE OF SOLAR MIN/MAX )
C        1988-11-26 ; T.KAJITA
C                     DATA ARE UP-DATED USING GAISSER AND VOLKOVA FLUX 
C                     WHICH INCLUDE MUON POLARIZATION EFFECT
C
C        1989-07-21 ; T.KAJITA
C                     LAST DATA UPDATE
C                     (A,B,BB,RLMIN,RLMAX,RMMIN,RMMAX)
C                      A AND B FOR NEU(E),NEU(MU),NEU(E)BAR,NEU(MU)BAR
C                      THERE ARE APPLOXIMATELY A FACTOR OF 2 MARGINE
C                      IN A AND BB.
C                      THEREFORE ONE NEED NOT CHANGE THE FOLLOWING
C                      THREE DATA WHEN ONE MADE A SMALL MODIFICATION.
C
C        1994-09-05 ; T.KAJITA
C                     E(MIN) FOR NUMU WAS CHANGED
C                      (FROM 0.12GEV TO 0.01GEV)
C
C        1995-01-07 ; T.KAJITA
C                     THE ONLY INPUT PARAMETERS ARE A AND B.
C                     ALMOST INDEPENDENT ON THE CHANGE OF FLUX AND
C                     CROSS SECTION.)
C                     (ALMOST = UPTO FACTOR OF TWO.)
C                     (FOR OLD RNNEUE SEE RNNEUEO1)
C
C        2007-11-10 ; T.TANAKA
C                     ADD UPMU MODE
C
C        2007-12-05 ; G.Mitsuka
C                     Maximum energy is expanded to 90TeV in FC/PC

#include "nesolact.h"
#ifdef UPMU_RCK
#include "upmu.h"
#endif



C      DIMENSION A(4),B(4)
       REAL A(4),B(4)
C      COMMON/MODGEN/MODGE1,MODGE2,THRTOT,ENEU,DIRNU(3)
C      COMMON /SOLACT/ ISOL, SOLACT
C
      DATA A/ 400., 600., 130., 250./
      DATA B/ -2.3, -2.0, -2.1, -1.9/


C
C
      IF(ISOL.EQ.0) SOLACT=0.5
      IF(SOLACT.LT.0.) SOLACT=0.
      IF(SOLACT.GT.1.) SOLACT=1.
C-------------------------------------------------------------------
C
C
C      IF(MODGE1.EQ.-1)GO TO 3000
C
      RNNEUE=0.
      I=0
      IF(IPAR.EQ.12)I=1
      IF(IPAR.EQ.14)I=2
      IF(IPAR.EQ.16)I=2
      IF(IPAR.EQ.-12)I=3
      IF(IPAR.EQ.-14)I=4
      IF(IPAR.EQ.-16)I=4
      IF(I.EQ.0)RETURN
C
C
C
#ifdef UPMU_RCK
      EMIN=UPMU_ETH
      EMAX=UPMU_MAXENE
#else

      EMIN=0.01
*     Supports up to 90TeV (2007,Dec. G.Mitsuka)
      EMAX=90000.
#endif

  100 F=RLU(DUM)
      ETEMP=F*(EMAX**(B(I)+1.)-EMIN**(B(I)+1.))+EMIN**(B(I)+1.)
      E=ETEMP**(1./(B(I)+1.))
      FEMAX=A(I)*(E**B(I))
      CALL NEFLUX(IPAR,E,FLXVAL)
      FE=FLXVAL*FNTOTPAU(IPAR,E)
C CHECK
      IF(FE.GT.FEMAX)
     &  PRINT *,' **** ERROR IN RNNEUE **** ENEU FEMAX FE ',E,FEMAX,FE
C
      RNFE=RLU(DUM)*FEMAX
      IF(RNFE.GT.FE)GO TO 100
      RNNEUE=E
C
C CHECK MODGE1
C
C      IF(MODGE1.EQ.1 .AND. RNNEUE.GT.5.) GO TO 100
C      IF(MODGE1.EQ.2 .AND. RNNEUE.LT.5.) GO TO 100
      RETURN
C
C ++ FIXED ENERGY
C
 3000 RNNEUE=ENEU/1000.
      RETURN
      END
C
C ++ BLOCK DATA
C
C      BLOCK DATA
C
C/MODGEN/MODGE1,MODGE2,THRTOT
C    MODGE1 : =0  GENERATE ALL NEUTRINO EVENTS
C             =1  GENERATE ONLY E<5GEV NEUTRINO EVENTS
C             =2  GENERATE ONLY E>5GEV NEUTRINO EVENTS
C    MODGE2 : =0  GENERATE CERENKOV FOR ALL EVENTS
C             =1  GENERATE CERENKOV ONLY FOR EXPECTED TOTAL P.E.<THRTOT
C    THRTOT :     LOOK MODGE2
C
C      COMMON/MODGEN/MODGE1,MODGE2,THRTOT,ENEU,DIRNU(3)
C      DATA MODGE1/1/
C      DATA MODGE2/0/
C      DATA THRTOT/10000./
C      END

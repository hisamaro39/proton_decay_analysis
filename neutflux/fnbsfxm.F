************************************************************************
*     ---------------------------------
      FUNCTION FNBSFXM(ENEU,DIR,IPAR,index)
*     ---------------------------------
C     
*     (Purpose)
C        Honda Differential Flux for mid energy
C
C     (Input)
C       ENEU   : Energy of Neutrino(in GeV)
C       ORGDIR : Direction of the neutrino
C       IPAR   : Flavor of neutrino
C       Index  :       
C                   1 : Flux = Energy^{-2.71}
C                   2 : Flux = Energy^{-2.74}
C
*     (Output)
C       Flux
C      
C     (Creation Date and Author)
C      1998.07     made by J.Kameda based on fnhdfxm.F
C                                             
************************************************************************
C
      Implicit NONE
#include "bessflx.h"      
#include "nesolact.h"
      
      REAL*4     FNBSFXM

      REAL*4     ENEU
      REAL*4     DIR(3)
      INTEGER*4  IPAR
      INTEGER*4  INDEX
      
      INTEGER*4  IEL,IEH
      INTEGER*4  ITL,ITH
      INTEGER*4  IPL,IPH

      REAL*4     RTL,RTH
      REAL*4     RPL,RPH

      REAL*4     FXLELTPL,FXLELTPH,FXHELTPL,FXHELTPH
      REAL*4     FXLEHTPL,FXLEHTPH,FXHEHTPL,FXHEHTPH
      
      REAL*4     FLXLELT,FLXLEHT
      REAL*4     FLXHELT,FLXHEHT
      REAL*4     FLXLE,FLXHE
      REAL*4     PHI,APHI

      REAL*4     A,B

      INTEGER*4  IPKIND
C      INTEGER*4  IS,IZ,IA,IE,ICOS


C---CHECK COMMON VAR(ISOL == 0 <--> SOLACT = 0.5)
      IF (ISOL.eq.0) SOLACT=0.50000

C---CHECK KIND
      IF (IPAR.eq.12)  IPKIND=1
      IF (IPAR.eq.-12) IPKIND=2
      IF (IPAR.eq.14)  IPKIND=3
      IF (IPAR.eq.-14) IPKIND=4

C---CHECK INDEX
      IF ((INDEX.ne.1).and.(INDEX.ne.2)) then
         write(*,*) "*** ERROR *** Stop executing ***"
         write(*,*) "***       in FNBSFXM.F       ***"
         write(*,*) "   Index = ",index
         write(*,*) "   Index should be  1 or 2      "
         STOP
      ENDIF

C---DECIDE ENERGY BIN
      IF (ENEU.lt.BESEM(1)) THEN
         write(*,*) "*** ERROR *** Stop executing ***"
         write(*,*) "***       in FNBSFXM.F       ***"
         write(*,*) "ENEU(",ENEU,") < ELOW(",BESEM(1),")"
         STOP
      ENDIF


      IF (ENEU.gt.BESEM(81)) THEN
         write(*,*) "*** ERROR *** Stop executing ***"
         write(*,*) "***       in FNBSFXM.F       ***"
         write(*,*) "ENEU(",ENEU,") > EHIGH(",BESEM(81),")"
         STOP
      ENDIF
      
      DO 5 IEH=1,81
         IF (ENEU.lt.BESEM(IEH)) GOTO 7
 5    CONTINUE
      IF (IEH.eq.82) THEN
C         write(*,*) "ERROR! IEH=82"
         IEH=81
      ENDIF
      
 7    IEL=IEH-1
      IF (IEH.eq.0) THEN
         write(*,*) "IEH=0"
         STOP
      ENDIF
         
C---DECIDE ANGLE BIN
      IF (abs(DIR(3)).gt.1) THEN
         write(*,*) "*** ERROR *** Stop executing ***"
         write(*,*) "***       in FNBSFXM.F       ***" 
         write(*,*) "abs(DIR(3))=abs(",DIR(3),") > 1."
         STOP
      ENDIF

      IF (DIR(3).le.-0.95) THEN
         ITL=1
         ITH=1
         GOTO 10
      ENDIF
      IF (DIR(3).gt.0.95) THEN
         ITL=20
         ITH=20
         GOTO 10
      ENDIF
      ITL=INT(DIR(3)*10.+10.5)
      ITH=ITL+1
      IF (ITH.eq.21) THEN 
         ITH=20
         GOTO 10
      ENDIF

 10   CONTINUE

      IF (DIR(1).EQ.0.0 .AND. DIR(2).EQ.0.0) THEN
         PHI = 0.0
      ELSE
         PHI = ATAN2(DIR(2),DIR(1))
      ENDIF
      
      IF (PHI.le.0.) THEN
         APHI = PHI + 2*3.141593
      ELSE
         APHI = PHI
      ENDIF

      IF (ABS(PHI).lt.3.141593/12.) THEN
         IPL=12
         IPH=1
         GOTO 15
      ENDIF

      IPL=INT( APHI/(3.141593/6.)-0.5 )+1
      IPH=IPL+1

C---CONSIDER EACH FLUX for LOW ENERGY AND HIGH ENERGY

 15   CONTINUE

      RTL     = REAL(ITL)*0.1-1.05
      RTH     = REAL(ITH)*0.1-1.05

      IF (IPL.eq.12) THEN
         RPL  = PHI-(-3.141593/12.)
         RPH  = (3.141593/12.)-PHI
      ELSE
         RPL  = APHI-(3.141593/6.*(REAL(IPL)-0.5))
         RPH  = (3.141593/6.*(REAL(IPH)-0.5))-APHI
      ENDIF

 17   FXLELTPL =  BFLXM(IEL,ITL,IPL,IPKIND,index)

      FXLELTPH =  BFLXM(IEL,ITL,IPH,IPKIND,index)

      FLXLELT  = (FXLELTPH-FXLELTPL)/(3.141593/6.)*RPL
     $          + FXLELTPL
      
      FXLEHTPL =  BFLXM(IEL,ITH,IPL,IPKIND,index)

      FXLEHTPH =  BFLXM(IEL,ITH,IPH,IPKIND,index)

      FLXLEHT  = (FXLEHTPH-FXLEHTPL)/(3.141593/6.)*RPL
     $          + FXLEHTPL
      
      IF (ITL.eq.ITH) THEN
         FLXLE=FLXLELT
         GOTO 20
      ENDIF

      FLXLE   = (FLXLEHT-FLXLELT)/(RTH-RTL)*(DIR(3)-RTL)+FLXLELT

 20   FXHELTPL =  BFLXM(IEH,ITL,IPL,IPKIND,index)

      FXHELTPH =  BFLXM(IEH,ITL,IPH,IPKIND,index)

      FLXHELT  = (FXHELTPH-FXHELTPL)/(3.141593/6.)*RPL
     $          + FXHELTPL
      
      FXHEHTPL =  BFLXM(IEH,ITH,IPL,IPKIND,index)

      FXHEHTPH =  BFLXM(IEH,ITH,IPH,IPKIND,index)

      FLXHEHT  = (FXHEHTPH-FXHEHTPL)/(3.141593/6.)*RPL
     $          + FXHEHTPL

      IF (ITL.eq.ITH) THEN
         FLXHE=FLXHELT
         GOTO 30
      ENDIF

      FLXHE   = (FLXHEHT-FLXHELT)/(RTH-RTL)*(DIR(3)-RTL)+FLXHELT


 30   B       = ALOG(FLXLE/FLXHE)/ALOG(BESEM(IEL)/BESEM(IEH))
      A       = FLXLE/BESEM(IEL)**B
      FNBSFXM = A*ENEU**B

      RETURN
      END

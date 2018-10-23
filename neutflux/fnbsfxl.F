************************************************************************
*     ---------------------------------
      FUNCTION FNBSFXL(ENEU,DIR,IPAR)
*     ---------------------------------
C     
*     (Purpose)
C        HONDA Differential Flux for low energy
C
C     (Input)
C       ENEU   : Energy of Neutrino(in GeV)
C       ORGDIR : Direction of the neutrino
C       IPAR   : Flavor of neutrino
C
*     (Output)
C       Flux (m^{-2}sr^{-1}sec^{-1}GeV^{-1})
C      
C     (Creation Date and Author)
C
C       1998.07 J.Kameda based on fnhdfxl.F
C                                             
************************************************************************

C
      Implicit NONE
#include "bessflx.h"      
#include "nesolact.h"
      
      REAL*4     FNBSFXL
      
      REAL*4     ENEU
      REAL*4     DIR(3)
      INTEGER*4  IPAR

      INTEGER*4  IEL,IEH
      INTEGER*4  ITL,ITH

      REAL*4     RTL,RTH

      REAL*4     FLXLELT,FLXLEHT
      REAL*4     FLXHELT,FLXHEHT
      REAL*4     FLXLE,FLXHE

      REAL*4     A,B

      INTEGER*4  IPKIND

C---CHECK COMMON VAR(ISOL == 0 <--> SOLACT = 0.5)
C---AT THIS MOMENT THIS VALUE is MEANINGLESS..
      IF (ISOL.eq.0) SOLACT=0.50000

C---CHECK KIND
      IF (IPAR.eq.12)  IPKIND=1
      IF (IPAR.eq.-12) IPKIND=2
      IF (IPAR.eq.14)  IPKIND=3
      IF (IPAR.eq.-14) IPKIND=4

C---DECIDE ENERGY BIN
      IF (ENEU.lt.BESEL(1)) THEN
         write(*,*) "*** ERROR *** Stop executing ***"
         write(*,*) "ENEU(",ENEU,") < ELOW(",BESEL(1),")"
         STOP
      ENDIF

      IF (ENEU.gt.BESEL(72)) THEN
         write(*,*) "*** ERROR *** Stop executing ***"
         write(*,*) "ENEU(",ENEU,") > EHIGH(",BESEL(72),")"
         STOP
      ENDIF
      
      DO 5 IEH=1,72
         IF (ENEU.lt.BESEL(IEH)) GOTO 7
 5    CONTINUE
 7    IEL=IEH-1
         
C---DECIDE ANGLE BIN
      IF (abs(DIR(3)).gt.1) THEN
         write(*,*) "*** ERROR *** Stop executing ***"
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

C---CONSIDER EACH FLUX for LOW ENERGY AND HIGH ENERGY

      RTL     = REAL(ITL)*0.1-1.05
      RTH     = REAL(ITH)*0.1-1.05

      FLXLELT = BFLXL(IEL,ITL,IPKIND,1)*(1.-SOLACT)
     $         +BFLXL(IEL,ITL,IPKIND,2)*SOLACT
      FLXLEHT = BFLXL(IEL,ITH,IPKIND,1)*(1.-SOLACT)
     $         +BFLXL(IEL,ITH,IPKIND,2)*SOLACT
C      FLXLELT = BFLXL(IEL,ITL,IPKIND)
C      FLXLEHT = BFLXL(IEL,ITH,IPKIND)


      IF (ITL.eq.ITH) THEN
         FLXLE=FLXLELT
         GOTO 20
      ENDIF

      FLXLE   = (FLXLEHT-FLXLELT)/(RTH-RTL)*(DIR(3)-RTL)+FLXLELT

 20   FLXHELT = BFLXL(IEH,ITL,IPKIND,1)*(1.-SOLACT)
     $         +BFLXL(IEH,ITL,IPKIND,2)*SOLACT
      FLXHEHT = BFLXL(IEH,ITH,IPKIND,1)*(1.-SOLACT)
     $         +BFLXL(IEH,ITH,IPKIND,2)*SOLACT
C 20   FLXHELT = BFLXL(IEH,ITL,IPKIND)
C      FLXHEHT = BFLXL(IEH,ITH,IPKIND)


      IF (ITL.eq.ITH) THEN
         FLXHE=FLXHELT
         GOTO 30
      ENDIF

      FLXHE   = (FLXHEHT-FLXHELT)/(RTH-RTL)*(DIR(3)-RTL)+FLXHELT

 30   B       = ALOG(FLXLE/FLXHE)/ALOG(BESEL(IEL)/BESEL(IEH))
      A       = FLXLE/BESEL(IEL)**B
      FNBSFXL = A*ENEU**B

      RETURN
      END

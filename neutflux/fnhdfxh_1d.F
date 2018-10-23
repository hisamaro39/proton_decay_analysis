************************************************************************
*     ---------------------------------
      FUNCTION FNHDFXH_1D(ENEU,DIR,IPAR)
*     ---------------------------------
C     
*     (Purpose)
C        Honda Differential Flux for High energy
C
C     (Input)
C       ENEU   : Energy of Neutrino(in GeV)
C       ORGDIR : Direction of the neutrino
C       IPAR   : Flavor of neutrino
C
*     (Output)
C       Flux
C      
C     (Creation Date and Author)
C       1996.??.?? ; Y.Hayato
C       1997.??.?? ; Y.Hayato ( debug : Direction is changed to the
C                                       detector coordinates)       
C                                             
************************************************************************
C
      Implicit NONE
#include "hondaflx_1d.h"      
#include "nesolact.h"
      
      REAL*4     FNHDFXH_1D
      
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
      IF (ENEU.lt.HONEH_1D(1)) THEN
         write(*,*) "*** ERROR *** Stop executing ***"
         write(*,*) "ENEU(",ENEU,") < ELOW(",HONEH_1D(1),")"
         STOP
      ENDIF
      
      IF (ENEU.ge.HONEH_1D(31)) THEN
         write(*,*) "*** ERROR *** Stop executing ***"
         write(*,*) "ENEU(",ENEU,") > EHIGH(",HONEH_1D(31),")"
         STOP
      ENDIF
      
      DO 5 IEH=1,31
         IF (ENEU.lt.HONEH_1D(IEH)) GOTO 7
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

C      FLXLELT = HFLXH(IEL,ITL,IPKIND,0)*(1.-SOLACT)
C     $         +HFLXH(IEL,ITL,IPKIND,1)*SOLACT
C      FLXLEHT = HFLXH(IEL,ITH,IPKIND,0)*(1.-SOLACT)
C     $         +HFLXH(IEL,ITH,IPKIND,1)*SOLACT
      FLXLELT = HFLXH_1D(IEL,ITL,IPKIND)
      FLXLEHT = HFLXH_1D(IEL,ITH,IPKIND)

      IF (ITL.eq.ITH) THEN
         FLXLE=FLXLELT
         GOTO 20
      ENDIF

      FLXLE   = (FLXLEHT-FLXLELT)/(RTH-RTL)*(DIR(3)-RTL)+FLXLELT

C 20   FLXHELT = HFLXH(IEH,ITL,IPKIND,0)*(1.-SOLACT)
C     $         +HFLXH(IEH,ITL,IPKIND,1)*SOLACT
c      FLXHEHT = HFLXH(IEH,ITH,IPKIND,0)*(1.-SOLACT)
C     $         +HFLXH(IEH,ITH,IPKIND,1)*SOLACT
 20   FLXHELT = HFLXH_1D(IEH,ITL,IPKIND)
      FLXHEHT = HFLXH_1D(IEH,ITH,IPKIND)

      IF (ITL.eq.ITH) THEN
         FLXHE=FLXHELT
         GOTO 30
      ENDIF

      FLXHE   = (FLXHEHT-FLXHELT)/(RTH-RTL)*(DIR(3)-RTL)+FLXHELT

 30   B       = ALOG(FLXLE/FLXHE)/ALOG(HONEH_1D(IEL)/HONEH_1D(IEH))
      A       = FLXLE/HONEH_1D(IEL)**B
      FNHDFXH_1D = A*ENEU**B

      RETURN
      END

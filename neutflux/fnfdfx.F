************************************************************************
*     ---------------------------------
      FUNCTION FNFDFX(ENEU,DIR,IPAR)
*     ---------------------------------
C     
*     (Purpose)
C        Fluka Differential Flux
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
C       2000.10.13 ; A.Kibayashi
C                                             
************************************************************************

C
      Implicit NONE
#include "flukaflx.h"      
#include "nesolact.h"
      
      REAL*4     FNFDFX
      
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
*--- fluka3D energy bin is 61 ---*
      IF (ENEU.lt.FONE(1)) THEN
         write(*,*) "*** ERROR *** Stop executing ***"
         write(*,*) "ENEU(",ENEU,") < ELOW(",FONE(1),")"
         STOP
      ENDIF

      IF (ENEU.gt.FONE(61)) THEN
         write(*,*) "*** ERROR *** Stop executing ***"
         write(*,*) "ENEU(",ENEU,") > EHIGH(",FONE(61),")"
         STOP
      ENDIF
      
      DO 5 IEH=1,61
         IF (ENEU.lt.FONE(IEH)) GOTO 7
 5    CONTINUE
 7    IEL=IEH-1

C---DECIDE ANGLE BIN
*--- fluka3D cosz bin is 40 ---*
      IF (abs(DIR(3)).gt.1) THEN
         write(*,*) "*** ERROR *** Stop executing ***"
         write(*,*) "abs(DIR(3))=abs(",DIR(3),") > 1."
         STOP
      ENDIF

      IF (DIR(3).le.-0.975) THEN
         ITL=1
         ITH=1
         GOTO 10
      ENDIF
      IF (DIR(3).gt.0.975) THEN
         ITL=40
         ITH=40
         GOTO 10
      ENDIF
      ITL=INT(DIR(3)*20.+20.5)  ! depends on no. of bins
      ITH=ITL+1
      IF (ITH.eq.41) THEN 
         ITH=40
         GOTO 10
      ENDIF

 10   CONTINUE

C---CONSIDER EACH FLUX for LOW ENERGY AND HIGH ENERGY

      RTL     = REAL(ITL)*0.05-1.025 ! depends on no. of bins
      RTH     = REAL(ITH)*0.05-1.025

*-ak- fluka01 with solar-min only
      FLXLELT = FFLX(IEL,ITL,IPKIND,1)*(1.-SOLACT)
     $         +FFLX(IEL,ITL,IPKIND,2)*SOLACT
      FLXLEHT = FFLX(IEL,ITH,IPKIND,1)*(1.-SOLACT)
     $         +FFLX(IEL,ITH,IPKIND,2)*SOLACT
C      FLXLELT = FFLX(IEL,ITL,IPKIND) ! flux low energy low theta
C      FLXLEHT = FFLX(IEL,ITH,IPKIND) ! flux low energy high theta

      IF (ITL.eq.ITH) THEN
         FLXLE=FLXLELT
         GOTO 20
      ENDIF

      FLXLE   = (FLXLEHT-FLXLELT)/(RTH-RTL)*(DIR(3)-RTL)+FLXLELT

 20   FLXHELT = FFLX(IEH,ITL,IPKIND,1)*(1.-SOLACT)
     $         +FFLX(IEH,ITL,IPKIND,2)*SOLACT
      FLXHEHT = FFLX(IEH,ITH,IPKIND,1)*(1.-SOLACT)
     $         +FFLX(IEH,ITH,IPKIND,2)*SOLACT
C 20   FLXHELT = FFLX(IEH,ITL,IPKIND) ! flux high energy low theta
C      FLXHEHT = FFLX(IEH,ITH,IPKIND) ! flux low energy high theta

      IF (ITL.eq.ITH) THEN
         FLXHE=FLXHELT
         GOTO 30
      ENDIF

      FLXHE   = (FLXHEHT-FLXHELT)/(RTH-RTL)*(DIR(3)-RTL)+FLXHELT

 30   B       = ALOG(FLXLE/FLXHE)/ALOG(FONE(IEL)/FONE(IEH))
      A       = FLXLE/FONE(IEL)**B
      FNFDFX = A*ENEU**B

      RETURN
      END

************************************************************************
*     ---------------------------------
      REAL FUNCTION FNGAIFX(ENEU,ORGDIR,IPAR)
*     ---------------------------------
C
*     (Purpose)
C       Honda Differential Flux for Low energy
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

      Implicit NONE
#include "gaisser.h"      
#include "nesolact.h"

c      REAL*4     FNGAIFX
      
      REAL*4     ENEU
      REAL*4     ORGDIR(3)
      INTEGER*4  IPAR,I

      INTEGER*4  IEL,IEH
      INTEGER*4  ITL,ITH

      REAL*4     RTL,RTH

      REAL*4     FLXLELT,FLXLEHT
      REAL*4     FLXHELT,FLXHEHT
      REAL*4     FLXLE,FLXHE

      REAL*4     A,B

      INTEGER*4  IPKIND

      REAL*4    DIR(3)

      INTEGER*4 ISTART
      DATA      ISTART/0/


C--- GEOMETRY CONVERSION ---
C--- Invert Z -> -Z
C--- rotate 40.583 degree in X-Y plane

      DIR(3)=-1.*ORGDIR(3)
      DIR(1)= COS(3.141593*2435./60./180.)*ORGDIR(1)
     $       +SIN(3.141593*2435./60./180.)*ORGDIR(2)
      DIR(2)=-SIN(3.141593*2435./60./180.)*ORGDIR(1)
     $       +COS(3.141593*2435./60./180.)*ORGDIR(2)
    
C---DATA has been stored?
      IF (ISTART.eq.0) THEN
         ISTART=1
         CALL NERDGAI
      ENDIF

C---CHECK COMMON VAR(ISOL == 0 <--> SOLACT = 0.5)
C---AT THIS MOMENT THIS VALUE is MEANINGLESS..
      IF (ISOL.eq.0) SOLACT=0.50000

C---CHECK KIND
      IF (IPAR.eq.12)  IPKIND=1
      IF (IPAR.eq.-12) IPKIND=2
      IF (IPAR.eq.14)  IPKIND=3
      IF (IPAR.eq.-14) IPKIND=4

C---DECIDE ENERGY BIN
      IF (ENEU.lt.EGAIS(1)) THEN
         write(*,*) "*** ERROR *** Stop executing ***"
         write(*,*) "ENEU(",ENEU,") < ELOW(",EGAIS(1),")"
         STOP
      ENDIF

      IF (ENEU.gt.EGAIS(63)) THEN
         write(*,*) "*** ERROR *** Stop executing ***"
         write(*,*) "ENEU(",ENEU,") > EHIGH(",EGAIS(63),")"
         STOP
      ENDIF
      
      DO 5 IEH=1,63
         IF (ENEU.lt.EGAIS(IEH)) GOTO 7
 5    CONTINUE
 7    IEL=IEH-1
         
C---DECIDE ANGLE BIN
      IF (abs(DIR(3)).gt.1) THEN
         write(*,*) "*** ERROR *** Stop executing ***"
         write(*,*) "abs(DIR(3))=abs(",DIR(3),") > 1."
         STOP
      ENDIF

      IF (DIR(3).eq.COSGAI(1)) THEN
         ITL=1
         ITH=1
         GOTO 11
      ENDIF
      IF (DIR(3).eq.COSGAI(12)) THEN
         ITL=12
         ITH=12
         GOTO 11
      ENDIF
      
      DO 10 I=1,12
         IF (DIR(3).ge.COSGAI(I)) THEN
            ITL=I
            ITH=I-1
            GOTO 11
         ENDIF
 10   CONTINUE

 11   CONTINUE
      IF (ITH.eq.0) THEN
         write(*,*) "ERROR!(ITH=0,DIR(3)=",DIR(3),")"
         STOP
      ENDIF

C---CONSIDER EACH FLUX for LOW ENERGY AND HIGH ENERGY

      RTL     = COSGAI(ITL)
      RTH     = COSGAI(ITH)

      FLXLELT = FLGAIS(1,IEL,ITL,IPKIND)*(1.-SOLACT)
     $         +FLGAIS(2,IEL,ITL,IPKIND)*SOLACT
      FLXLEHT = FLGAIS(1,IEL,ITH,IPKIND)*(1.-SOLACT)
     $         +FLGAIS(2,IEL,ITH,IPKIND)*SOLACT

      IF (ITL.eq.ITH) THEN
         FLXLE=FLXLELT
         GOTO 20
      ENDIF

      FLXLE   = (FLXLEHT-FLXLELT)/(RTH-RTL)*(DIR(3)-RTL)+FLXLELT

 20   FLXHELT = FLGAIS(1,IEH,ITL,IPKIND)*(1.-SOLACT)
     $         +FLGAIS(2,IEH,ITL,IPKIND)*SOLACT
      FLXHEHT = FLGAIS(1,IEH,ITH,IPKIND)*(1.-SOLACT)
     $         +FLGAIS(2,IEH,ITH,IPKIND)*SOLACT

      IF (ITL.eq.ITH) THEN
         FLXHE=FLXHELT
         GOTO 30
      ENDIF

      FLXHE   = (FLXHEHT-FLXHELT)/(RTH-RTL)*(DIR(3)-RTL)+FLXHELT

 30   B       = ALOG(FLXLE/FLXHE)/ALOG(EGAIS(IEL)/EGAIS(IEH))
      A       = FLXLE/EGAIS(IEL)**B
      FNGAIFX = A*ENEU**B

      RETURN
      END

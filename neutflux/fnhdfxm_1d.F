************************************************************************
*     ---------------------------------
      FUNCTION FNHDFXM_1D(ENEU,DIR,IPAR)
*     ---------------------------------
C     
*     (Purpose)
C        Honda Differential Flux for mid energy
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
      
      REAL*4     FNHDFXM_1D

      REAL*4     ENEU
      REAL*4     DIR(3)
      INTEGER*4  IPAR

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

C---DECIDE ENERGY BIN
      IF (ENEU.lt.HONEM_1D(1)) THEN
         write(*,*) "*** ERROR *** Stop executing ***"
         write(*,*) "ENEU(",ENEU,") < ELOW(",HONEM_1D(1),")"
         STOP
      ENDIF

      IF (ENEU.gt.HONEM_1D(50)) THEN
         write(*,*) "*** ERROR *** Stop executing ***"
         write(*,*) "ENEU(",ENEU,") > EHIGH(",HONEM_1D(50),")"
         STOP
      ENDIF
      
      DO 5 IEH=1,50
         IF (ENEU.lt.HONEM_1D(IEH)) GOTO 7
 5    CONTINUE
      IF (IEH.eq.51) THEN
C         write(*,*) "ERROR! IEH=51"
         IEH=50
      ENDIF
      
 7    IEL=IEH-1
      IF (IEH.eq.0) THEN
         write(*,*) "IEH=0"
         STOP
      ENDIF
         
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

      
C      DO 900 IS=1,2
C         DO 920 IZ=1,20
C            DO 930 IA=1,12
C               DO 910 IE=1,50
C                  IF (IZ.le.11) THEN
C                     ICOS=IZ+9
C                  ELSE
C                     ICOS=21-IZ
C                  ENDIF
C                  write(*,'(E11.4,E11.4,E11.4,E11.4,E11.4)')
C     $                 HONEM(IE),
C     $                 HFLXM(IE,ICOS,IA,3,IS),
C     $                 HFLXM(IE,ICOS,IA,4,IS),
C     $                 HFLXM(IE,ICOS,IA,1,IS),
C     $                 HFLXM(IE,ICOS,IA,2,IS)
C 910           CONTINUE
C               
C 930        CONTINUE
C 920     CONTINUE
C 900  CONTINUE

C      write(*,*) "HONEM(",IEL,")=",HONEM(IEL)," / ENEU=",ENEU,
C     $           " / HONEM(",IEH,")=",HONEM(IEH)
      
C      write(*,*) "HFLXM(",IEL,",",ITL,",",IPL,",",IPKIND,",",1,")=",
C     $     HFLXM(IEL,ITL,IPL,IPKIND,1)
C      write(*,*) "HFLXM(",IEL,",",ITH,",",IPL,",",IPKIND,",",1,")=",
C     $     HFLXM(IEL,ITH,IPL,IPKIND,1)

 17   FXLELTPL =( HFLXM_1D(IEL,ITL,IPL,IPKIND,1)*(1.-SOLACT)
     $           +HFLXM_1D(IEL,ITL,IPL,IPKIND,2)*SOLACT    )
      FXLELTPH =( HFLXM_1D(IEL,ITL,IPH,IPKIND,1)*(1.-SOLACT)
     $           +HFLXM_1D(IEL,ITL,IPH,IPKIND,2)*SOLACT    )
      FLXLELT  = (FXLELTPH-FXLELTPL)/(3.141593/6.)*RPL
     $          + FXLELTPL
      
      FXLEHTPL =( HFLXM_1D(IEL,ITH,IPL,IPKIND,1)*(1.-SOLACT)
     $           +HFLXM_1D(IEL,ITH,IPL,IPKIND,2)*SOLACT    )
      FXLEHTPH =( HFLXM_1D(IEL,ITH,IPH,IPKIND,1)*(1.-SOLACT)
     $           +HFLXM_1D(IEL,ITH,IPH,IPKIND,2)*SOLACT    )
      FLXLEHT  = (FXLEHTPH-FXLEHTPL)/(3.141593/6.)*RPL
     $          + FXLEHTPL
      
      IF (ITL.eq.ITH) THEN
         FLXLE=FLXLELT
         GOTO 20
      ENDIF

      FLXLE   = (FLXLEHT-FLXLELT)/(RTH-RTL)*(DIR(3)-RTL)+FLXLELT

 20   FXHELTPL =( HFLXM_1D(IEH,ITL,IPL,IPKIND,1)*(1.-SOLACT)
     $           +HFLXM_1D(IEH,ITL,IPL,IPKIND,2)*SOLACT    )
      FXHELTPH =( HFLXM_1D(IEH,ITL,IPH,IPKIND,1)*(1.-SOLACT)
     $           +HFLXM_1D(IEH,ITL,IPH,IPKIND,2)*SOLACT    )
      FLXHELT  = (FXHELTPH-FXHELTPL)/(3.141593/6.)*RPL
     $          + FXHELTPL
      
      FXHEHTPL =( HFLXM_1D(IEH,ITH,IPL,IPKIND,1)*(1.-SOLACT)
     $           +HFLXM_1D(IEH,ITH,IPL,IPKIND,2)*SOLACT    )
      FXHEHTPH =( HFLXM_1D(IEH,ITH,IPH,IPKIND,1)*(1.-SOLACT)
     $           +HFLXM_1D(IEH,ITH,IPH,IPKIND,2)*SOLACT    )
      FLXHEHT  = (FXHEHTPH-FXHEHTPL)/(3.141593/6.)*RPL
     $          + FXHEHTPL

      IF (ITL.eq.ITH) THEN
         FLXHE=FLXHELT
         GOTO 30
      ENDIF

      FLXHE   = (FLXHEHT-FLXHELT)/(RTH-RTL)*(DIR(3)-RTL)+FLXHELT

 30   B       = ALOG(FLXLE/FLXHE)/ALOG(HONEM_1D(IEL)/HONEM_1D(IEH))
      A       = FLXLE/HONEM_1D(IEL)**B
      FNHDFXM_1D = A*ENEU**B

      RETURN
      END

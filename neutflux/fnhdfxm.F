************************************************************************
*     ---------------------------------
      FUNCTION FNHDFXM(ENEU,DIR,IPAR)
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
C       2001.09.18 ; M.Ishitsuka 
C           replace flux tables 
C                  honda97mid.dat & honda96high.dat -> hkkm2.dat
C                                             
C       2007.12.05 ; G.Mitsuka, FC/PC's binning is merged with upmu
C
C       2008.01.21 ; G.Mitsuka, Re-weighting Volkova flux above 10TeV
C
************************************************************************
C
      Implicit NONE
#include "hondaflx.h"      
#include "nesolact.h"
#include "hondaver.h"
#include "necardap.h"
      
      REAL*4     FNHDFXM

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

      REAL*4 HondaR
C      INTEGER*4  IS,IZ,IA,IE,ICOS


C---CHECK COMMON VAR(ISOL == 0 <--> SOLACT = 0.5)
      IF (ISOL.eq.0) SOLACT=0.50000

C---CHECK KIND
      IF (IPAR.eq.12)  IPKIND=1
      IF (IPAR.eq.-12) IPKIND=2
      IF (IPAR.eq.14)  IPKIND=3
      IF (IPAR.eq.-14) IPKIND=4

C---DECIDE ENERGY BIN
      IF (ENEU.lt.HONEM(1)) THEN
         write(*,*) "*** ERROR *** Stop executing ***"
         write(*,*) "ENEU(",ENEU,") < ELOW(",HONEM(1),")"
         STOP
      ENDIF

      IF (ENEU.gt.HONEM(107)) THEN
         write(*,*) "*** ERROR *** Stop executing ***"
         write(*,*) "ENEU(",ENEU,") > EHIGH(",HONEM(107),")"
         STOP
      ENDIF
      
      DO 5 IEH=1,107
         IF (ENEU.lt.HONEM(IEH)) GOTO 7
 5    CONTINUE
      IF (IEH.eq.108) THEN
         IEH=107
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

C     ------------------------------------------------------------
C     Absolute value of Volkova flux is adjusted 
C     following the version of Honda flux.
C     Then if this routine is called by other version's fnhonfx??,
C     Volkova flux must be re-adjusted.

      if (ENEU .gt. 1.0E+04) then ! Ratio is 1. below 10TeV

         if ( (NEFLXFLG.eq.0 .and. HONDAVER.eq.1 ) .or.
     $        (NEFLXFLG.ne.0 .and. HONDAVER.eq.0 ) ) then
C     Common value and reading code are same version
            HondaR = 1.
         elseif ( NEFLXFLG.eq.0 .and. HONDAVER.eq.0 ) then
C     e.g. Common value is Honda06(NEFLXFLG=0) and HONDAVER is Honda03(HONDAVER=0)
            if (IPKIND.eq.1) HondaR = 1.6E-2/0.03039
            if (IPKIND.eq.2) HondaR = 2.1E-2/0.03375
            if (IPKIND.eq.3) HondaR = 5.4E-1/0.96184
            if (IPKIND.eq.4) HondaR = 5.7E-1/0.84404
            
         elseif ( NEFLXFLG.ne.0 .and. HONDAVER.eq.1 ) then
C     e.g. Common value is Honda03(NEFLXFLG=6) and HONDAVER is Honda06(HONDAVER=1)
            if (IPKIND.eq.1) HondaR = 0.03039/1.6E-2
            if (IPKIND.eq.2) HondaR = 0.03375/2.1E-2
            if (IPKIND.eq.3) HondaR = 0.96184/5.4E-1
            if (IPKIND.eq.4) HondaR = 0.84404/5.7E-1
            
         else
            print*, 'NEFLXFLG is ',NEFLXFLG,' HONDAVER is ',HONDAVER
            print*, 'This mode is not supportted'
            stop
         endif
      else
         HondaR = 1.
      endif
C     ------------------------------------------------------------

      
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

 17   FXLELTPL =( HFLXM(IEL,ITL,IPL,IPKIND,1)*(1.-SOLACT)
     $           +HFLXM(IEL,ITL,IPL,IPKIND,2)*SOLACT    )*HondaR
      FXLELTPH =( HFLXM(IEL,ITL,IPH,IPKIND,1)*(1.-SOLACT)
     $           +HFLXM(IEL,ITL,IPH,IPKIND,2)*SOLACT    )*HondaR
      FLXLELT  = (FXLELTPH-FXLELTPL)/(3.141593/6.)*RPL
     $          + FXLELTPL
      
      FXLEHTPL =( HFLXM(IEL,ITH,IPL,IPKIND,1)*(1.-SOLACT)
     $           +HFLXM(IEL,ITH,IPL,IPKIND,2)*SOLACT    )*HondaR
      FXLEHTPH =( HFLXM(IEL,ITH,IPH,IPKIND,1)*(1.-SOLACT)
     $           +HFLXM(IEL,ITH,IPH,IPKIND,2)*SOLACT    )*HondaR
      FLXLEHT  = (FXLEHTPH-FXLEHTPL)/(3.141593/6.)*RPL
     $          + FXLEHTPL
      
      IF (ITL.eq.ITH) THEN
         FLXLE=FLXLELT
         GOTO 20
      ENDIF

      FLXLE   = (FLXLEHT-FLXLELT)/(RTH-RTL)*(DIR(3)-RTL)+FLXLELT

 20   FXHELTPL =( HFLXM(IEH,ITL,IPL,IPKIND,1)*(1.-SOLACT)
     $           +HFLXM(IEH,ITL,IPL,IPKIND,2)*SOLACT    )*HondaR
      FXHELTPH =( HFLXM(IEH,ITL,IPH,IPKIND,1)*(1.-SOLACT)
     $           +HFLXM(IEH,ITL,IPH,IPKIND,2)*SOLACT    )*HondaR
      FLXHELT  = (FXHELTPH-FXHELTPL)/(3.141593/6.)*RPL
     $          + FXHELTPL
      
      FXHEHTPL =( HFLXM(IEH,ITH,IPL,IPKIND,1)*(1.-SOLACT)
     $           +HFLXM(IEH,ITH,IPL,IPKIND,2)*SOLACT    )*HondaR
      FXHEHTPH =( HFLXM(IEH,ITH,IPH,IPKIND,1)*(1.-SOLACT)
     $           +HFLXM(IEH,ITH,IPH,IPKIND,2)*SOLACT    )*HondaR
      FLXHEHT  = (FXHEHTPH-FXHEHTPL)/(3.141593/6.)*RPL
     $          + FXHEHTPL

      IF (ITL.eq.ITH) THEN
         FLXHE=FLXHELT
         GOTO 30
      ENDIF

      FLXHE   = (FLXHEHT-FLXHELT)/(RTH-RTL)*(DIR(3)-RTL)+FLXHELT

 30   B       = ALOG(FLXLE/FLXHE)/ALOG(HONEM(IEL)/HONEM(IEH))
      A       = FLXLE/HONEM(IEL)**B
      FNHDFXM = A*ENEU**B

      RETURN
      END

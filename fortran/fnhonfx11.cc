      REAL FUNCTION FNHONFX11(ENUE,ORGDIR,SOLACT,IPAR)

c Return differential Honda flux GeV-1/m^2/s/sr
c
c Input:
c     ENUE    Energy of neutrino (GeV)
c     ORGDIR  Direction of neutrino origin (Super-K convention)
c             orgdir(3)<0 : downward direction
c             orgdir(3)>0 : upward direction
c     SOLACT  Solar activity (0. = min, 1. = max)
c     IPAR    Flavor of neutrino(12=nue,-12=nuebar,14=numu,-14=numubar)
c
c     Nov 2011    Tristan Irvine
c                 based on fnhonfx06.F (G. Mitsuka)
c
c     flux tables: kam-solmax-mountain.d; kam-solmin-mountain.d
c		available at http://www.icrr.u-tokyo.ac.jp/~mhonda/nflx2011/index.html

      IMPLICIT NONE
#include "hondaver.h"
c
c     Indices are: Azimuth(12),  Zenith(20),  Energy(101),  Solar Min/Max(2)
c     ranges are : 15 - 345, 0.95 - -0.95, 1.0E-01 - 1.0E+04, Solar min/max
c

      INTEGER*4 LUNHFLX
      PARAMETER (LUNHFLX=68)

      
      REAL ENERGY(101)
      REAL BFNU(101,20,12,4,2)
c ----------------------------------------------------------------------
      REAL FINT
      INTEGER*4 IHNDL
      EXTERNAL FINT
      REAL FNHONFX
      EXTERNAL FNHONFX
      REAL FNHONFX11LOW
      EXTERNAL FNHONFX11LOW
      CHARACTER*132 DUMMY_STRING, FILENAME
      REAL ENUE, ORGDIR(3), DIR(3)
      INTEGER IPAR 
      INTEGER J, NLINES
      REAL BIN, SOLACT
      INTEGER IAZI, ICOSZEN, ILOGE, ISOLACT
      LOGICAL INITIALIZED
      DATA INITIALIZED/.FALSE./
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
      print *, 'Start fnhonfx11 !!'
      print *, "ENUE=", ENUE
      print *, "SOLACT=", SOLACT
      print *, "IPAR=", IPAR
      print *, "ORGDIR 1/2/3=",ORGDIR(1),"/",ORGDIR(2),"/",ORGDIR(3)
c ======================================================================
c
c     Initialize upon first call by loading tables:
c     =============================================
      
	  IF (.NOT.INITIALIZED) THEN
         WRITE(*,'(A)') 
     &        " Now LOADING HONDA-2011 FLUX"
c
c        Flux tables file (solar min & solar max together):
c        --------------------
         CALL SKOPENF(LUNHFLX,1,'f',IHNDL)
         IF (IHNDL.LT.0) THEN
            print *, 'Could not open hkkm11mt.dat flux file'
            stop
         ENDIF
c     ------------------solact-----------------------------------
         DO ISOLACT=1,2
c     ------------------zenith-----------------------------------
            DO ICOSZEN = 1, 20
c     ------------------azimth-----------------------------------
               DO IAZI = 1, 12
                  READ(LUNHFLX,'(A)') DUMMY_STRING
                  
                  READ(LUNHFLX,'(A)') DUMMY_STRING
                  

c     ------------------energy-----------------------------------
                  DO ILOGE = 1, 101
                     READ(LUNHFLX,*) ENERGY(ILOGE),
     &                    BFNU(ILOGE,21-ICOSZEN,IAZI,3,ISOLACT),
     &                    BFNU(ILOGE,21-ICOSZEN,IAZI,4,ISOLACT),
     &                    BFNU(ILOGE,21-ICOSZEN,IAZI,1,ISOLACT),
     &                    BFNU(ILOGE,21-ICOSZEN,IAZI,2,ISOLACT)
                  END DO
               END DO
            END DO
         END DO
         CALL SKCLOSEF(LUNHFLX)
         INITIALIZED = .TRUE.
      END IF                    ! End of first-time initialization
C---CHECK KIND
      IF (IPAR.eq.12)  IPKIND=1
      IF (IPAR.eq.-12) IPKIND=2
      IF (IPAR.eq.14)  IPKIND=3
      IF (IPAR.eq.-14) IPKIND=4

c
c     Geometry conversion, made to agree with neut/fnhonfx.F
c     ================================================
      DIR(3) = -ORGDIR(3)
      DIR(1)=-(COS(3.141593*2435./60./180.)*ORGDIR(1)
     &       +SIN(3.141593*2435./60./180.)*ORGDIR(2))
      DIR(2)=-(-SIN(3.141593*2435./60./180.)*ORGDIR(1)
     &       +COS(3.141593*2435./60./180.)*ORGDIR(2))
c
c     At high energy or low energy, different routines are required:
c     ================================
      IF ((ENUE.GT.ENERGY(101)).or.(ENUE.LT.ENERGY(1))) THEN
         HONDAVER = 1

         IF ((ENUE.LT.ENERGY(1))) THEN
c     If energy is <0.1GeV, use the 2011 Honda Low-e flux
             FNHONFX11 = FNHONFX11LOW(ENUE,ORGDIR,SOLACT,IPAR)
         ELSE 
c     If energy is >10TeV, use the 2001 Honda flux
             FNHONFX11 = FNHONFX(ENUE,ORGDIR,IPAR)  
         ENDIF
         HONDAVER = 0           ! re-set as default
      ELSE IF ((ENUE.LE.ENERGY(101)) .AND. (ENUE.GE.ENERGY(1))) THEN
c     Here the correct bins are being allocated for the function's input values.
         DO 5 IEH=1,101
            IF (ENUE.lt.ENERGY(IEH)) GOTO 7
 5       CONTINUE
c     The following if statement is included to stop the code breaking at exactly 10TeV.
 7       IF ((IEH.gt.101)) THEN
             IEH=101
         ENDIF 
         IEL=IEH-1
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
         
 10      CONTINUE
         
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

 15      CONTINUE

         RTL     = REAL(ITL)*0.1-1.05
         RTH     = REAL(ITH)*0.1-1.05
         
         IF (IPL.eq.12) THEN
            RPL  = PHI-(-3.141593/12.)
            RPH  = (3.141593/12.)-PHI
         ELSE
            RPL  = APHI-(3.141593/6.*(REAL(IPL)-0.5))
            RPH  = (3.141593/6.*(REAL(IPH)-0.5))-APHI
         ENDIF
c Below section is extrapolating between theta/thi/energy bins,
c and solar min/maximum, to get intermediate values of flux.         
 17      FXLELTPL =( BFNU(IEL,ITL,IPL,IPKIND,1)*(1.-SOLACT)
     $        +BFNU(IEL,ITL,IPL,IPKIND,2)*SOLACT    )
         FXLELTPH =( BFNU(IEL,ITL,IPH,IPKIND,1)*(1.-SOLACT)
     $        +BFNU(IEL,ITL,IPH,IPKIND,2)*SOLACT    )
         FLXLELT  = (FXLELTPH-FXLELTPL)/(3.141593/6.)*RPL
     $        + FXLELTPL
         
         FXLEHTPL =( BFNU(IEL,ITH,IPL,IPKIND,1)*(1.-SOLACT)
     $        +BFNU(IEL,ITH,IPL,IPKIND,2)*SOLACT    )
         FXLEHTPH =( BFNU(IEL,ITH,IPH,IPKIND,1)*(1.-SOLACT)
     $        +BFNU(IEL,ITH,IPH,IPKIND,2)*SOLACT    )
         FLXLEHT  = (FXLEHTPH-FXLEHTPL)/(3.141593/6.)*RPL
     $        + FXLEHTPL
         
         IF (ITL.eq.ITH) THEN
            FLXLE=FLXLELT
            GOTO 20
         ENDIF
         
         FLXLE   = (FLXLEHT-FLXLELT)/(RTH-RTL)*(DIR(3)-RTL)+FLXLELT
 20      FXHELTPL =( BFNU(IEH,ITL,IPL,IPKIND,1)*(1.-SOLACT)
     $        +BFNU(IEH,ITL,IPL,IPKIND,2)*SOLACT    )
         FXHELTPH =( BFNU(IEH,ITL,IPH,IPKIND,1)*(1.-SOLACT)
     $        +BFNU(IEH,ITL,IPH,IPKIND,2)*SOLACT    )
         FLXHELT  = (FXHELTPH-FXHELTPL)/(3.141593/6.)*RPL
     $          + FXHELTPL
         
         FXHEHTPL =( BFNU(IEH,ITH,IPL,IPKIND,1)*(1.-SOLACT)
     $        +BFNU(IEH,ITH,IPL,IPKIND,2)*SOLACT    )
         FXHEHTPH =( BFNU(IEH,ITH,IPH,IPKIND,1)*(1.-SOLACT)
     $        +BFNU(IEH,ITH,IPH,IPKIND,2)*SOLACT    )
         FLXHEHT  = (FXHEHTPH-FXHEHTPL)/(3.141593/6.)*RPL
     $        + FXHEHTPL
         
         IF (ITL.eq.ITH) THEN
            FLXHE=FLXHELT
            GOTO 30
         ENDIF
         
         FLXHE   = (FLXHEHT-FLXHELT)/(RTH-RTL)*(DIR(3)-RTL)+FLXHELT
         
 30      B       = ALOG(FLXLE/FLXHE)/ALOG(ENERGY(IEL)/ENERGY(IEH))
         A       = FLXLE/ENERGY(IEL)**B
         FNHONFX11 = A*ENUE**B

      END IF
      print *,'FNHONFX11=' ,FNHONFX11

      RETURN
      END

      

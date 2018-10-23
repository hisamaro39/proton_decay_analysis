      REAL FUNCTION FNFLUFX03(ENUE,ORGDIR,SOLACT,IPAR)

c Return differential Fluka flux GeV-1/m^2/s/sr
c
c Input:
c     ENEU    Energy of neutrino (GeV)
c     ORGDIR  Direction of neutrino origin (Super-K convention)
c     SOLACT  Solar activity (0. = min, 1. = max)
c     IPAR    Flavor of neutrino(12=nue,-12=nuebar,14=numu,-14=numubar)
c
c     Mar 2004    Masaki Ishitsuka 
c                 based on fnhonfx.F (Y.Hayato)
c
c     format of the table : fluka03.dat
c ============================================================================
c  unit    /GeV^1/m^2/s/sr
c ============================================================================
c  comment (smim or smax)
c  comment (cos=*)
c  comment (azimth=*)
c  number  (energy nue sig_nue nuebar sig_nuebar num sig_num numbar sig_numbar)
c ============================================================================
      IMPLICIT NONE
c
c     Indices are: Azimuth(6),  Zenith(20),  Energy(100),  Solar Min/Max(2)
c     ranges are : 30 - 330, -0.95 - 0.95, 0.106 - 9440.61, Min - Max
c
      INTEGER*4 LUNFFLX
      PARAMETER (LUNFFLX=63)

      INTEGER*4 IHNDL

      REAL LOGE(100),ENERGY(100)
      REAL BFNU(100,20,6,4,2)
c ----------------------------------------------------------------------
      REAL FINT
      EXTERNAL FINT
      REAL FNHONFX
      EXTERNAL FNHONFX
      CHARACTER*132 DUMMY_STRING, FILENAME
      REAL ENUE, ORGDIR(3), DIR(3)
      INTEGER IPAR 
      INTEGER J
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
c ======================================================================
c
c     Initialize upon first call by loading tables:
c     =============================================
      IF (.NOT.INITIALIZED) THEN
         WRITE(*,'(A)') 
     &        " Now LOADING FLUKA-2003 FLUX"
c
c        Load solar min file:
c        --------------------
         CALL SKOPENF(LUNFFLX,1,'f',IHNDL)
         IF (IHNDL.LT.0) THEN
            print *, 'Could not open fluka03 flux file'
            stop
         ENDIF
c         CALL FINDCONSTS(FILENAME,'Fluka Flux 2003',
c     &        'fluka03.dat',0)
c         OPEN(63,FILE=FILENAME,FORM='FORMATTED',STATUS='OLD')
c     ------------------solact-----------------------------------
         DO ISOLACT=1,2
            READ(LUNFFLX,'(A)') DUMMY_STRING
c     ------------------zenith-----------------------------------
            DO ICOSZEN = 1, 20
               READ(LUNFFLX,'(A)') DUMMY_STRING
c     ------------------azimth-----------------------------------
               DO IAZI = 1, 6
                  READ(LUNFFLX,'(A)') DUMMY_STRING
c     ------------------energy-----------------------------------
                  DO ILOGE = 1, 100
                     READ(LUNFFLX,*) LOGE(ILOGE),
     &                    BFNU(ILOGE,ICOSZEN,IAZI,1,ISOLACT),
     &                    BFNU(ILOGE,ICOSZEN,IAZI,2,ISOLACT),
     &                    BFNU(ILOGE,ICOSZEN,IAZI,3,ISOLACT),
     &                    BFNU(ILOGE,ICOSZEN,IAZI,4,ISOLACT)
                     ENERGY(ILOGE)=10**LOGE(ILOGE)
                  END DO
               END DO
            END DO
         END DO
         CALL SKCLOSEF(LUNFFLX)
c         CLOSE(63)
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
c     At high energy or low energy, use honda-3d routine:
c     ================================
      IF ((ENUE.GT.10**LOGE(100)).or.(ENUE.LT.10**LOGE(1))) THEN
         FNFLUFX03 = FNHONFX(ENUE,ORGDIR,IPAR)
      ELSE IF ((ENUE.LE.10**LOGE(100)) .AND. (ENUE.GE.10**LOGE(1))) THEN
         DO 5 IEH=1,100
            IF (ENUE.lt.ENERGY(IEH)) GOTO 7
 5       CONTINUE
 7       IEL=IEH-1

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
         
         IF (ABS(PHI).lt.3.141593/6.) THEN
            IPL=6
            IPH=1
            GOTO 15
         ENDIF

         IPL=INT( APHI/(3.141593/3.)-0.5 )+1
         IPH=IPL+1

C---CONSIDER EACH FLUX for LOW ENERGY AND HIGH ENERGY

 15      CONTINUE

         RTL     = REAL(ITL)*0.1-1.05
         RTH     = REAL(ITH)*0.1-1.05
         
         IF (IPL.eq.12) THEN
            RPL  = PHI-(-3.141593/6.)
            RPH  = (3.141593/6.)-PHI
         ELSE
            RPL  = APHI-(3.141593/3.*(REAL(IPL)-0.5))
            RPH  = (3.141593/3.*(REAL(IPH)-0.5))-APHI
         ENDIF
         
 17      FXLELTPL =( BFNU(IEL,ITL,IPL,IPKIND,1)*(1.-SOLACT)
     $        +BFNU(IEL,ITL,IPL,IPKIND,2)*SOLACT    )
         FXLELTPH =( BFNU(IEL,ITL,IPH,IPKIND,1)*(1.-SOLACT)
     $        +BFNU(IEL,ITL,IPH,IPKIND,2)*SOLACT    )
         FLXLELT  = (FXLELTPH-FXLELTPL)/(3.141593/3.)*RPL
     $        + FXLELTPL
         
         FXLEHTPL =( BFNU(IEL,ITH,IPL,IPKIND,1)*(1.-SOLACT)
     $        +BFNU(IEL,ITH,IPL,IPKIND,2)*SOLACT    )
         FXLEHTPH =( BFNU(IEL,ITH,IPH,IPKIND,1)*(1.-SOLACT)
     $        +BFNU(IEL,ITH,IPH,IPKIND,2)*SOLACT    )
         FLXLEHT  = (FXLEHTPH-FXLEHTPL)/(3.141593/3.)*RPL
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
         FLXHELT  = (FXHELTPH-FXHELTPL)/(3.141593/3.)*RPL
     $          + FXHELTPL
         
         FXHEHTPL =( BFNU(IEH,ITH,IPL,IPKIND,1)*(1.-SOLACT)
     $        +BFNU(IEH,ITH,IPL,IPKIND,2)*SOLACT    )
         FXHEHTPH =( BFNU(IEH,ITH,IPH,IPKIND,1)*(1.-SOLACT)
     $        +BFNU(IEH,ITH,IPH,IPKIND,2)*SOLACT    )
         FLXHEHT  = (FXHEHTPH-FXHEHTPL)/(3.141593/3.)*RPL
     $        + FXHEHTPL
         
         IF (ITL.eq.ITH) THEN
            FLXHE=FLXHELT
            GOTO 30
         ENDIF
         FLXHE   = (FLXHEHT-FLXHELT)/(RTH-RTL)*(DIR(3)-RTL)+FLXHELT
         
 30      B       = ALOG(FLXLE/FLXHE)/ALOG(ENERGY(IEL)/ENERGY(IEH))
         A       = FLXLE/ENERGY(IEL)**B
         FNFLUFX03 = A*ENUE**B

      END IF

      RETURN
      END

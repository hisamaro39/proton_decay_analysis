      REAL FUNCTION FNBARFX03(ENUE,ORGDIR,SOLACT,IPAR)

c Return differential Bartol flux
c
c Input:
c     ENEU    Energy of neutrino (GeV)
c     ORGDIR  Direction of neutrino origin (Super-K convention)
c     SOLACT  Solar activity (0. = min, 1. = max)
c     IPAR    Flavor of neutrino(12=nue,-12=nuebar,14=numu,-14=numubar)
c
c Masaki Ishitsuka, Nov 2003 but the table is up to 10 GeV.
c                            It will be updated soon.
c
c     2004 .01.29; C.Saji
c     The table is updated up to 8892.0112 GeV and format is changed.
c     2004. 03.19; M.Ishitsuka
c     Replace flux table to updated one: bartol03.dat
c     2012. 01.23; K.Okumura
c     bug fixed to avoid accessing out-of-bound array
c
c     format of the table : bartol03.dat
c ============================================================================
c  unit    (E*dN/dE) /cm^2/s/sr
c ============================================================================
c  comment (smim or smax)
c  comment (cos=*)
c  comment (azimth=*)
c  number  (energy nue sig_nue nuebar sig_nuebar num sig_num numbar sig_numbar)
c ============================================================================
      IMPLICIT NONE
c
c     Indices are: Azimuth(18), Zenith(20), Energy(40), Solar Min/Max(2)
c     ranges are : 10 - 350, -0.95 - 0.95, 0.106 - 9.441, Min - Max
c
      INTEGER*4 LUNBFLX,NPARA
      PARAMETER (LUNBFLX=62)
      PARAMETER (NPARA=70)

      INTEGER*4 IHNDL

      REAL LOGE(NPARA),ENERGY(NPARA)
      REAL BFNU(NPARA,20,18,4,2)
      REAL SIGNU(NPARA,20,18,4,2)
c ----------------------------------------------------------------------
      REAL FINT
      EXTERNAL FINT
      REAL FNHONFX
      EXTERNAL FNHONFX
      CHARACTER*132 DUMMY_STRING
      REAL DUMMY_REAL(9)
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
     &        " Now LOADING BARTOL-2003 FLUX"
c
c        Load solar min file:
c        --------------------
         CALL SKOPENF(LUNBFLX,1,'f',IHNDL)
         IF (IHNDL.LT.0) THEN
            print *, 'Could not open bartol03 flux file'
            stop
         ENDIF
c         CALL FINDCONSTS(FILENAME,'Bartol Flux 2003',
c     &        'bartol03.dat',0)
c         OPEN(62,FILE=FILENAME,FORM='FORMATTED',STATUS='OLD')
c     ------------------solact-----------------------------------
         DO ISOLACT=1,2
            READ(LUNBFLX,'(A)') DUMMY_STRING
c     ------------------zenith-----------------------------------
            DO ICOSZEN = 1, 20
               READ(LUNBFLX,'(A)') DUMMY_STRING
c     ------------------azimth-----------------------------------
               DO IAZI = 1, 18
                  READ(LUNBFLX,'(A)') DUMMY_STRING
c     ------------------energy-----------------------------------
                  DO ILOGE = 1, 70
                     IF(ILOGE.GT.NPARA) THEN
                        READ(LUNBFLX,*) DUMMY_REAL
                     ELSE
                        READ(LUNBFLX,*) LOGE(ILOGE),
     &                       BFNU(ILOGE,ICOSZEN,IAZI,1,ISOLACT),
     &                       SIGNU(ILOGE,ICOSZEN,IAZI,1,ISOLACT),
     &                       BFNU(ILOGE,ICOSZEN,IAZI,2,ISOLACT),
     &                       SIGNU(ILOGE,ICOSZEN,IAZI,2,ISOLACT),
     &                       BFNU(ILOGE,ICOSZEN,IAZI,3,ISOLACT),
     &                       SIGNU(ILOGE,ICOSZEN,IAZI,3,ISOLACT),
     &                       BFNU(ILOGE,ICOSZEN,IAZI,4,ISOLACT),
     &                       SIGNU(ILOGE,ICOSZEN,IAZI,4,ISOLACT)
                        ENERGY(ILOGE)=10**LOGE(ILOGE)
                     ENDIF
                  END DO
               END DO
            END DO
         END DO
         CALL SKCLOSEF(LUNBFLX)
c         CLOSE(62)
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
      IF ((ENUE.GT.10**LOGE(NPARA)).or.(ENUE.LT.10**LOGE(1))) THEN
         FNBARFX03 = FNHONFX(ENUE,ORGDIR,IPAR)
      ELSE IF ((ENUE.LE.10**LOGE(NPARA))
     $        .AND. (ENUE.GE.10**LOGE(1))) THEN
         DO 5 IEH=1,70
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

CKO     commint out          
CKO         IF (ABS(PHI).lt.3.141593/18.) THEN
CKO            IPL=18
CKO            IPH=1
CKO            GOTO 15
CKO         ENDIF

         IPL=INT( APHI/(3.141593/9.)-0.5 )+1
         IPH=IPL+1

CKO     added for avoiding out of bounds
         IF (ABS(PHI).lt.3.141593/18.) THEN
            IPL=18
            IPH=1
         ENDIF
         IF ( IPL.EQ.0 ) THEN
            IPL=18
         ENDIF
         IF ( IPL.EQ.18 ) THEN
            IPH=1
         ENDIF

         if ( ipl.gt.18 ) then
            ipl=18
            iph=1
         endif
         if ( iph.gt.18 ) then
            ipl=17
            iph=18
         endif


C---CONSIDER EACH FLUX for LOW ENERGY AND HIGH ENERGY

 15      CONTINUE

         RTL     = REAL(ITL)*0.1-1.05
         RTH     = REAL(ITH)*0.1-1.05

CKO     fixed bug         
CKO         IF (IPL.eq.12) THEN
CKO            RPL  = PHI-(-3.141593/18.)
CKO            RPH  = (3.141593/12.)-PHI
         IF (IPL.eq.18) THEN
            RPL  = PHI-(-3.141593/18.)
            RPH  = (3.141593/18.)-PHI
         ELSE
            RPL  = APHI-(3.141593/9.*(REAL(IPL)-0.5))
            RPH  = (3.141593/9.*(REAL(IPH)-0.5))-APHI
         ENDIF

C         write(99,*)ipl,iph,phi,aphi

         
 17      FXLELTPL =( BFNU(IEL,ITL,IPL,IPKIND,1)*(1.-SOLACT)
     $        +BFNU(IEL,ITL,IPL,IPKIND,2)*SOLACT    )
         FXLELTPH =( BFNU(IEL,ITL,IPH,IPKIND,1)*(1.-SOLACT)
     $        +BFNU(IEL,ITL,IPH,IPKIND,2)*SOLACT    )
         FLXLELT  = (FXLELTPH-FXLELTPL)/(3.141593/9.)*RPL
     $        + FXLELTPL
         
         FXLEHTPL =( BFNU(IEL,ITH,IPL,IPKIND,1)*(1.-SOLACT)
     $        +BFNU(IEL,ITH,IPL,IPKIND,2)*SOLACT    )
         FXLEHTPH =( BFNU(IEL,ITH,IPH,IPKIND,1)*(1.-SOLACT)
     $        +BFNU(IEL,ITH,IPH,IPKIND,2)*SOLACT    )
         FLXLEHT  = (FXLEHTPH-FXLEHTPL)/(3.141593/9.)*RPL
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
         FLXHELT  = (FXHELTPH-FXHELTPL)/(3.141593/9.)*RPL
     $          + FXHELTPL
         
         FXHEHTPL =( BFNU(IEH,ITH,IPL,IPKIND,1)*(1.-SOLACT)
     $        +BFNU(IEH,ITH,IPL,IPKIND,2)*SOLACT    )
         FXHEHTPH =( BFNU(IEH,ITH,IPH,IPKIND,1)*(1.-SOLACT)
     $        +BFNU(IEH,ITH,IPH,IPKIND,2)*SOLACT    )
         FLXHEHT  = (FXHEHTPH-FXHEHTPL)/(3.141593/9.)*RPL
     $        + FXHEHTPL
         
         IF (ITL.eq.ITH) THEN
            FLXHE=FLXHELT
            GOTO 30
         ENDIF
         
         FLXHE   = (FLXHEHT-FLXHELT)/(RTH-RTL)*(DIR(3)-RTL)+FLXHELT
         
 30      B       = ALOG(FLXLE/FLXHE)/ALOG(ENERGY(IEL)/ENERGY(IEH))
         A       = FLXLE/ENERGY(IEL)**B
         FNBARFX03 = A*ENUE**B/ENUE

      END IF

      RETURN
      END

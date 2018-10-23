      REAL FUNCTION FNBARFX(ENUE,ORGDIR,SOLACT,IPAR)
c ======================================================================
c Return differential Bartol flux (E*dN/dE) /cm^2/s/sr
c
c Input:
c     ENEU    Energy of neutrino (GeV)
c     ORGDIR  Direction of neutrino origin (Super-K convention)
c     SOLACT  Solar activity (0. = min, 1. = max)
c     IPAR    Flavor of neutrino(12=nue,-12=nuebar,14=numu,-14=numubar)
c
c Ed Kearns, April 1998
c ======================================================================
      IMPLICIT NONE
c
c     Indices are: Azimuth(14), Zenith(20), Energy(23), Solar Min/Max(2)
c     Note: Azimuth has 2 extra bins to aid interpolation at 0/360.
c
      REAL LOGE(23)
      REAL BFNUE(14,20,23,2), BFNUEBAR(14,20,23,2)
      REAL BFNUMU(14,20,23,2), BFNUMUBAR(14,20,23,2)
c ----------------------------------------------------------------------
      REAL FINT
      EXTERNAL FINT
      REAL FNGAIFX
      EXTERNAL FNGAIFX
      CHARACTER*132 DUMMY_STRING, FILENAME
      REAL ENUE, ORGDIR(3), DIR(3)
      INTEGER IPAR 
      INTEGER J, NA(4)
      REAL BIN, SOLACT, XINT(4), A(59)
      INTEGER IAZI, ICOSZEN, ILOGE, ISOLACT
      LOGICAL INITIALIZED
      DATA INITIALIZED/.FALSE./
      DATA NA/14,20,23,2/
c ======================================================================
c
c     Initialize upon first call by loading tables:
c     =============================================
      IF (.NOT.INITIALIZED) THEN
c
c        Initialize array of grid-points for the interpolation routine:
c        --------------------------------------------------------------
         J = 1
c
c        Add 1 extra bin, before and after, in phi, 
c        to wrap-around interpolation at phi = 0/360.
c        So note that NA(1) is 14 and not 12.
c        ------------------------------------
         BIN = 2.*3.1415926536/12.
         A(J) = -.5*BIN
         J = J + 1
         DO IAZI = 1, 12
            A(J) = FLOAT(IAZI-1)*BIN + .5*BIN
            J = J + 1
         END DO
         A(J) = A(J-1) + BIN
         J = J + 1
c
         BIN = .1
         DO ICOSZEN = 1, 20
            A(J) = -0.95 + FLOAT(ICOSZEN-1)*BIN
            J = J + 1
         END DO
c
         BIN = .1
         DO ILOGE = 1, 23
            A(J) = -1.25 + FLOAT(ILOGE-1)*BIN ! GeV
            J = J + 1
         END DO
c
         BIN = 1.
         DO ISOLACT = 1, 2
            A(J) = FLOAT(ISOLACT - 1)*BIN
            J = J + 1
         END DO

         WRITE(*,'(A)') 
     &        " Now LOADING PHI-DEPENDENT BARTOL FLUX"
c
c        Load solar min file:
c        --------------------
         CALL FINDCONSTS(FILENAME,'Bartol Flux Solar Min',
     &        'kam_phi.smin',0)
         OPEN(61,FILE=FILENAME,FORM='FORMATTED',STATUS='OLD')
         DO ICOSZEN = 1, 10
            READ(61,'(A)') DUMMY_STRING
c 
c           Save IAZI=1,14 to wrap around interpolation at 0/360.
c           -----------------------------------------------------
            DO IAZI = 2, 13
               READ(61,'(A)') DUMMY_STRING
 
               DO ILOGE = 1, 23
c
c                 Bartol flux tables have one line with both +cosz and -cosz
c                 ----------------------------------------------------------
                  READ(61,*) LOGE(ILOGE),
     &                    BFNUE(IAZI,21-ICOSZEN,ILOGE,1),
     &                    BFNUEBAR(IAZI,21-ICOSZEN,ILOGE,1),
     &                    BFNUMU(IAZI,21-ICOSZEN,ILOGE,1),
     &                    BFNUMUBAR(IAZI,21-ICOSZEN,ILOGE,1),
     &                    BFNUE(IAZI,ICOSZEN,ILOGE,1),
     &                    BFNUEBAR(IAZI,ICOSZEN,ILOGE,1),
     &                    BFNUMU(IAZI,ICOSZEN,ILOGE,1),
     &                    BFNUMUBAR(IAZI,ICOSZEN,ILOGE,1)
               END DO
            END DO
         END DO
c
c        Fill in azimuth bins that interpolate at 0/360 degrees:
c        -------------------------------------------------------
         DO ICOSZEN = 1, 20
            DO ILOGE = 1, 23
               BFNUE(1,ICOSZEN,ILOGE,1) = 
     &              BFNUE(13,ICOSZEN,ILOGE,1)
               BFNUEBAR(1,ICOSZEN,ILOGE,1) = 
     &              BFNUEBAR(13,ICOSZEN,ILOGE,1)
               BFNUMU(1,ICOSZEN,ILOGE,1) = 
     &              BFNUMU(13,ICOSZEN,ILOGE,1)
               BFNUMUBAR(1,ICOSZEN,ILOGE,1) =
     &              BFNUMUBAR(13,ICOSZEN,ILOGE,1)
               BFNUE(14,ICOSZEN,ILOGE,1) = 
     &              BFNUE(2,ICOSZEN,ILOGE,1)
               BFNUEBAR(14,ICOSZEN,ILOGE,1) = 
     &              BFNUEBAR(2,ICOSZEN,ILOGE,1)
               BFNUMU(14,ICOSZEN,ILOGE,1) = 
     &              BFNUMU(2,ICOSZEN,ILOGE,1)
               BFNUMUBAR(14,ICOSZEN,ILOGE,1) =
     &              BFNUMUBAR(2,ICOSZEN,ILOGE,1)
            END DO
         END DO
         CLOSE(61)
c
c        Load solar max file:
c        --------------------
         CALL FINDCONSTS(FILENAME,'Bartol Flux Solar Max',
     &        'kam_phi.smax',0)
         OPEN(61,FILE=FILENAME,FORM='FORMATTED',STATUS='OLD')
         DO ICOSZEN = 1, 10
            READ(61,'(A)') DUMMY_STRING
c 
c           Save IAZI=1,14 to wrap around interpolation at 0/360.
c           -----------------------------------------------------
            DO IAZI = 2, 13
               READ(61,'(A)') DUMMY_STRING
 
               DO ILOGE = 1, 23
c
c                 Bartol flux tables have one line with both +cosz and -cosz
c                 ----------------------------------------------------------
                  READ(61,*) LOGE(ILOGE),
     &                    BFNUE(IAZI,21-ICOSZEN,ILOGE,2),
     &                    BFNUEBAR(IAZI,21-ICOSZEN,ILOGE,2),
     &                    BFNUMU(IAZI,21-ICOSZEN,ILOGE,2),
     &                    BFNUMUBAR(IAZI,21-ICOSZEN,ILOGE,2),
     &                    BFNUE(IAZI,ICOSZEN,ILOGE,2),
     &                    BFNUEBAR(IAZI,ICOSZEN,ILOGE,2),
     &                    BFNUMU(IAZI,ICOSZEN,ILOGE,2),
     &                    BFNUMUBAR(IAZI,ICOSZEN,ILOGE,2)
               END DO
            END DO
         END DO
c
c        Fill in azimuth bins that interpolate at 0/360 degrees:
c        -------------------------------------------------------
         DO ICOSZEN = 1, 20
            DO ILOGE = 1, 23
               BFNUE(1,ICOSZEN,ILOGE,2) = 
     &              BFNUE(13,ICOSZEN,ILOGE,2)
               BFNUEBAR(1,ICOSZEN,ILOGE,2) = 
     &              BFNUEBAR(13,ICOSZEN,ILOGE,2)
               BFNUMU(1,ICOSZEN,ILOGE,2) = 
     &              BFNUMU(13,ICOSZEN,ILOGE,2)
               BFNUMUBAR(1,ICOSZEN,ILOGE,2) =
     &              BFNUMUBAR(13,ICOSZEN,ILOGE,2)
               BFNUE(14,ICOSZEN,ILOGE,2) = 
     &              BFNUE(2,ICOSZEN,ILOGE,2)
               BFNUEBAR(14,ICOSZEN,ILOGE,2) = 
     &              BFNUEBAR(2,ICOSZEN,ILOGE,2)
               BFNUMU(14,ICOSZEN,ILOGE,2) = 
     &              BFNUMU(2,ICOSZEN,ILOGE,2)
               BFNUMUBAR(14,ICOSZEN,ILOGE,2) =
     &              BFNUMUBAR(2,ICOSZEN,ILOGE,2)
            END DO
         END DO
         CLOSE(61)
         INITIALIZED = .TRUE.
      END IF                    ! End of first-time initialization
c
c     Geometry conversion, made to agree with neut/fnhonfx.F
c     Note: a minus sign was introduced between 97f and 98a,
c     this convention agrees with 98a (and later). 
c     Some MC files were generated with 97f however...
c     ================================================
      DIR(3) = -ORGDIR(3)
      DIR(1)=(COS(3.141593*2435./60./180.)*ORGDIR(1)
     &       +SIN(3.141593*2435./60./180.)*ORGDIR(2))
      DIR(2)=(-SIN(3.141593*2435./60./180.)*ORGDIR(1)
     &       +COS(3.141593*2435./60./180.)*ORGDIR(2))
c
c     At high energy, use old routine:
c     ================================
      IF (ENUE .GT. 10**LOGE(23)) THEN
         FNBARFX = FNGAIFX(ENUE,DIR,IPAR)
      ELSE IF ((ENUE.LE.10**LOGE(23)) .AND. (ENUE.GE.10**LOGE(1))) THEN
c
c        Multi-dimensional linear interpolation (CERNLIB):
c        =================================================
         XINT(1) = ATAN2(DIR(2),DIR(1))
         IF (XINT(1) .LT. 0) XINT(1) = 2.*3.1415926536 + XINT(1)
         XINT(2) = DIR(3)
         XINT(3) = LOG10(ENUE)
         XINT(4) = SOLACT
         IF (IPAR .EQ. 12) THEN
            FNBARFX = FINT(4,XINT,NA,A,BFNUE)
         ELSE IF (IPAR .EQ. -12) THEN
            FNBARFX = FINT(4,XINT,NA,A,BFNUEBAR)
         ELSE IF (IPAR .EQ. 14) THEN
            FNBARFX = FINT(4,XINT,NA,A,BFNUMU)
         ELSE IF (IPAR .EQ. -14) THEN
            FNBARFX = FINT(4,XINT,NA,A,BFNUMUBAR)
         ELSE
            WRITE(*,'(A,1X,I3)') 
     &           "FNBARFX=0.0 called with unknown particle type",IPAR
            FNBARFX = 0.0
         END IF
c
c        Convert from E * cm-2 s-1 sr-1 to m-2 s-1 sr-1
c        to match convention of fnhonfx and fngaifx:
c        ===========================================
         FNBARFX = FNBARFX*10000./ENUE
      ELSE
         WRITE(*,'(A,1X,G9.4)') 
     &        "FNBARFX=0.0 called with too low energy =",ENUE
         FNBARFX = 0.0
      END IF

      RETURN
      END

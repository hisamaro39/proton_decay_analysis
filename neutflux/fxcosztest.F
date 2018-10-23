      PROGRAM FXCOSZTEST
c ----------------------------------------------------------------------
c 
c     Load Bartol and Honda flux tables and check zenith flux.
c     Output is written to fxcosztest.out.
c     Results may be plotted with fxcosztest.kumac
c
c     Ed Kearns, April 1998
c ----------------------------------------------------------------------
      IMPLICIT NONE
      INCLUDE 'nesolact.h'
      REAL FNBARFX, FNHONFX
      EXTERNAL FNBARFX, FNHONFX
      REAL DEGTORAD
      PARAMETER (DEGTORAD = 3.1415926536/180.)
      REAL E, COSZEN, PHI, LOGE, DIR(3)
      REAL BMUFXSMIN, BMUFXSMAX, BEFXSMIN, BEFXSMAX
      REAL HMUFXSMIN, HMUFXSMAX, HEFXSMIN, HEFXSMAX
      INTEGER IAZI, ICOSZEN, ILOGE
c ----------------------------------------------------------------------
      OPEN(33,FILE='fxcosztest.out',FORM='FORMATTED',STATUS='UNKNOWN')
      DO ILOGE = 5, 23, 2
         LOGE = -1.25 + (ILOGE-1)*.1
         E = 10**LOGE
         WRITE(33,317) E
         DO ICOSZEN = 20, 1, -1
            COSZEN = 1.05 - .1*FLOAT(ICOSZEN)
            BMUFXSMIN = 0.0
            BMUFXSMAX = 0.0
            BEFXSMIN  = 0.0
            BEFXSMAX  = 0.0
            HMUFXSMIN = 0.0
            HMUFXSMAX = 0.0
            HEFXSMIN  = 0.0
            HEFXSMAX  = 0.0
            DO IAZI = 1, 12
               PHI = 15. + FLOAT(IAZI-1)*30.
               DIR(1) = -COS(PHI*DEGTORAD)
               DIR(2) = -SIN(PHI*DEGTORAD)
               DIR(3) = -COSZEN
               BMUFXSMIN = BMUFXSMIN + 
     &                     FNBARFX(E,DIR,0.,14) +
     &                     FNBARFX(E,DIR,0.,-14)
               BMUFXSMAX = BMUFXSMAX + 
     &                     FNBARFX(E,DIR,1.,14) +
     &                     FNBARFX(E,DIR,1.,-14)
               BEFXSMIN = BEFXSMIN + 
     &                    FNBARFX(E,DIR,0.,12) +
     &                    FNBARFX(E,DIR,0.,-12)
               BEFXSMAX = BEFXSMAX + 
     &                    FNBARFX(E,DIR,1.,12) +
     &                    FNBARFX(E,DIR,1.,-12)

               ISOL = 1
               SOLACT = 0.   ! set solar activity via common
               HMUFXSMIN = HMUFXSMIN + 
     &                     FNHONFX(E,DIR,14) +
     &                     FNHONFX(E,DIR,-14)
               SOLACT = 1.
               HMUFXSMAX = HMUFXSMAX + 
     &                     FNHONFX(E,DIR,14) +
     &                     FNHONFX(E,DIR,-14)
               SOLACT = 0.
               HEFXSMIN = HEFXSMIN + 
     &                    FNHONFX(E,DIR,12) +
     &                    FNHONFX(E,DIR,-12)
               SOLACT = 1.
               HEFXSMAX = HEFXSMAX + 
     &                    FNHONFX(E,DIR,12) +
     &                    FNHONFX(E,DIR,-12)
            END DO
c
c           Change units to E*dN/dE cm-2 s-1 sr-1
c
            BMUFXSMIN = BMUFXSMIN*E/10000.
            BMUFXSMAX = BMUFXSMAX*E/10000.
            BEFXSMIN = BEFXSMIN*E/10000.
            BEFXSMAX = BEFXSMAX*E/10000.

            HMUFXSMIN = HMUFXSMIN*E/10000.
            HMUFXSMAX = HMUFXSMAX*E/10000.
            HEFXSMIN = HEFXSMIN*E/10000.
            HEFXSMAX = HEFXSMAX*E/10000.

            WRITE(33,318) COSZEN,
     &           BMUFXSMIN/12.,BMUFXSMAX/12.,BEFXSMIN/12.,BEFXSMAX/12.,
     &           HMUFXSMIN/12.,HMUFXSMAX/12.,HEFXSMIN/12.,HEFXSMAX/12.
         END DO
      END DO
c
c     Now do some high energy points
c
      DO ILOGE = 0, 8, 2
         LOGE = 1.15 + FLOAT(ILOGE)*.1
         E = 10**LOGE
         WRITE(33,317) E
         DO ICOSZEN = 20, 1, -1
            COSZEN = 1.05 - .1*FLOAT(ICOSZEN)
            BMUFXSMIN = 0.0
            BMUFXSMAX = 0.0
            BEFXSMIN  = 0.0
            BEFXSMAX  = 0.0
            HMUFXSMIN = 0.0
            HMUFXSMAX = 0.0
            HEFXSMIN  = 0.0
            HEFXSMAX  = 0.0
            DO IAZI = 1, 12
               PHI = 15. + FLOAT(IAZI-1)*30.
               DIR(1) = -COS(PHI*DEGTORAD)
               DIR(2) = -SIN(PHI*DEGTORAD)
               DIR(3) = -COSZEN
               BMUFXSMIN = BMUFXSMIN + 
     &                     FNBARFX(E,DIR,0.,14) +
     &                     FNBARFX(E,DIR,0.,-14)
               BMUFXSMAX = BMUFXSMAX + 
     &                     FNBARFX(E,DIR,1.,14) +
     &                     FNBARFX(E,DIR,1.,-14)
               BEFXSMIN = BEFXSMIN + 
     &                    FNBARFX(E,DIR,0.,12) +
     &                    FNBARFX(E,DIR,0.,-12)
               BEFXSMAX = BEFXSMAX + 
     &                    FNBARFX(E,DIR,1.,12) +
     &                    FNBARFX(E,DIR,1.,-12)

               ISOL = 1
               SOLACT = 0.   ! set solar activity via common
               HMUFXSMIN = HMUFXSMIN + 
     &                     FNHONFX(E,DIR,14) +
     &                     FNHONFX(E,DIR,-14)
               SOLACT = 1.
               HMUFXSMAX = HMUFXSMAX + 
     &                     FNHONFX(E,DIR,14) +
     &                     FNHONFX(E,DIR,-14)
               SOLACT = 0.
               HEFXSMIN = HEFXSMIN + 
     &                    FNHONFX(E,DIR,12) +
     &                    FNHONFX(E,DIR,-12)
               SOLACT = 1.
               HEFXSMAX = HEFXSMAX + 
     &                    FNHONFX(E,DIR,12) +
     &                    FNHONFX(E,DIR,-12)
            END DO
c
c           Change units to E*dN/dE cm-2 s-1 sr-1
c
            BMUFXSMIN = BMUFXSMIN*E/10000.
            BMUFXSMAX = BMUFXSMAX*E/10000.
            BEFXSMIN = BEFXSMIN*E/10000.
            BEFXSMAX = BEFXSMAX*E/10000.

            HMUFXSMIN = HMUFXSMIN*E/10000.
            HMUFXSMAX = HMUFXSMAX*E/10000.
            HEFXSMIN = HEFXSMIN*E/10000.
            HEFXSMAX = HEFXSMAX*E/10000.

            WRITE(33,318) COSZEN,
     &           BMUFXSMIN/12.,BMUFXSMAX/12.,BEFXSMIN/12.,BEFXSMAX/12.,
     &           HMUFXSMIN/12.,HMUFXSMAX/12.,HEFXSMIN/12.,HEFXSMAX/12.
         END DO
      END DO
317   FORMAT('*** Neutrino energy = ',F8.3,' GeV')
318   FORMAT(1X,F5.2,8(1X,G9.3))
      END

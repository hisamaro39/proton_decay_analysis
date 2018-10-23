*******************************************************************
*      -----------------------------------------
       SUBROUTINE NELEEFLX(IPAR,ENEU,DIR,FLUXNE)
*      -----------------------------------------
*
*      ( Purpose )
*        GET THE ATMOSPHERIC NEUTRINO FLUX FOR A GIVEN NENU,DIR,IPAR
*        THE FLUX IS THE ONE CALCULATED BY LEE IN 1989
*
*      ( Input )
*        IPAR; NEUTRINO TYPE
*             12; NUE
*            -12; NUEBAR
*             14; NUMU
*            -14; NUMUBAR
*        ENEU; NEUTRINO ENERGY IN GEV
*        DIR ; NEUTRINO DIRECTION ( IF DIR=(0,0,0) FLUXNE = Integrated flux )
*
*      ( Output )
*        FLUXNE; NEUTRINO FLUX     ( /M**2/SEC/SR/GEV )
*
*     ( Creation Date and Auther )
*       1988.11.24 ; T.KAJITA
*                    LATEST HONDA-FLUX IS INSTALLED (<7GEV)  (1990)
*                    VOLKOVA FLUX IS INSTALLED FOR E(NEU)>10GEV
*                    NEU/(NEU-BAR) RATIO IS INSTALLED ACCORDING TO THE 
*                    HONDA'S TABLE
*                    ALSO NEU/(NEU-BAR) AT HIGH ENERGY IS ESTIMATED 
*                    FROM THE SAME TABLE
*                    CORRECTION OF NEU(E)/NEU(MU) RATIO AT LOW E (<100MEV) 
*                    IS ESTIMATED FROM THE HIGHER ENERGY CALCULATION
*       1996.01.24 ; Y.KANAYA
*                    FOR SK
*
*      ( Comment )
*        FOR THE ANGULAR DEPENDENCE OF THE FLUX, WE USED THAT OF HONDA.
*        SOLAR MIN,MAX OR ANY COMBINATION OF MIN/MAX IS NOT AVAILABLE
*            IF(SOLACT=1.) -- SOLAR MAX.                  =MEAN
*            IF(SOLACT=0.) -- SOLAR MIN.                  =MEAN
*         (  IF(SOLACT=0.5) -- AVERAGE OF SOLAR MIN/MAX ) =MEAN
*        CALL NERDLEE & NERDHON BEFOREHAND
*
***********************************************************************
#include "nesolact.h"
      DIMENSION DIRCAL(3), DIR(3), ENEHI(3), PCOSH(8),
     &         FEHI(8,3),  FMHI(8,3)
      COMMON / NELEEF / ENELOW(135),
     &       FEMIN(135), FEMAX(135), FEBMIN(135), FEBMAX(135),
     &       FMMIN(135), FMMAX(135), FMBMIN(135), FMBMAX(135)
C
      DATA ENEHI  / 10., 30., 100. /
      DATA PCOSH  /1.0, 0.60, 0.40, 0.30, 0.20, 0.10, 0.05, 0.00 /
C
C NEU(E)+NEU(E)-BAR FLUX AT HIGH ENERGIES (VOLKOVA) /M**2/SEC/SR/GEV
      DATA FEHI/ .681E-1,  1.17E-1,  1.69E-1,  2.13E-1,  2.83E-1,
     &           3.66E-1,  3.99E-1,  4.18E-1,
     &           1.45E-3,  2.49E-3,  3.74E-3,  4.90E-3,  7.00E-3,
     &           10.1E-3,  11.6E-3,  12.7E-3,
     &           2.15E-5,  3.40E-5,  5.00E-5,  6.70E-5,  10.1E-5,
     &           16.5E-5,  20.6E-5,  24.2E-5          /
C
C NEU(MU)+NEU(MU)-BAR FLUX AT HIGH ENERGIES (VOLKOVA) /M**2/SEC/SR/GEV
      DATA FMHI/ 4.38E-1,  5.19E-1,  5.91E-1,  6.50E-1,  7.33E-1,
     &           8.28E-1,  8.65E-1,  8.86E-1,
     &           1.53E-2,  1.86E-2,  2.17E-2,  2.42E-2,  2.79E-2,
     &           3.25E-2,  3.44E-2,  3.57E-2,
     &           3.91E-4,  4.95E-4,  5.89E-4,  6.62E-4,  7.70E-4,
     &           9.19E-4,  9.91E-4,  10.5E-4           /
C
C -- STATEMENT FUNCTION
C             NEU/(NEU+NEUBAR) AT E>7GEV FOR NEU(E) AND NEU(MU)
      FNENEB(COSTH)=0.5575
      FNMNMB(COSTH)=0.5025+ 0.045*ABS(COSTH)
      DATA ICALL/0/
C
C --  READ TABLE
C
      IF (ICALL.EQ.0) THEN
         ICALL = 1
         CALL NERDLEE
      ENDIF
C
C
      IF(ISOL.EQ.0) SOLACT=0.5
      IF(SOLACT.LT.0.) SOLACT=0.
      IF(SOLACT.GT.1.) SOLACT=1.
C
C
C --DETERMINE ENERGY BINS
C
      DO 1 I=1,125
      IMAX=I
      DENEU = ENEU-ENELOW(I)
      IF(DENEU.LT.0.) THEN
           ILOW=I-1
           IF(ILOW.EQ.0) ILOW=1
           IHIGH=I
           IENERG=1
           GO TO 10
           END IF
    1 CONTINUE
C
      DO 2 J=1,3
CCCCC DO 2 J=2,3
      JMAX=J
      DENEU = ENEU-ENEHI(J)
      IF(DENEU.LT.0.) THEN
           ILOW=J-1
           IHIGH=J
           IENERG=3
           IF(ILOW.EQ.0) THEN
CCCCC      IF(ILOW.EQ.1) THEN
                         ILOW=IMAX
                         IENERG=2
                         END IF
           GO TO 10
           END IF
    2 CONTINUE
           ILOW=JMAX
           IHIGH=JMAX
           IENERG=3
   10 CONTINUE
C
C
      GO TO (100,200,300 ),IENERG
C
C
C -- E(NEU) LESS THAN 3 GEV
  100 CONTINUE
C
      IF(IPAR.EQ. 12) THEN
             FLCOSL= SOLACT*( FEMAX( ILOW ) +FEMAX(ILOW ))/2. +
     &          (1.-SOLACT)*( FEMIN( ILOW ) +FEMIN(ILOW ))/2.
             FLCOSH= SOLACT*( FEMAX( IHIGH) +FEMAX(IHIGH))/2. +
     &          (1.-SOLACT)*( FEMIN( IHIGH) +FEMIN(IHIGH))/2.
             END IF
      IF(IPAR.EQ.-12) THEN
             FLCOSL= SOLACT*( FEBMAX( ILOW ) +FEBMAX(ILOW ))/2. +
     &          (1.-SOLACT)*( FEBMIN( ILOW ) +FEBMIN(ILOW ))/2.
             FLCOSH= SOLACT*( FEBMAX( IHIGH) +FEBMAX(IHIGH))/2. +
     &          (1.-SOLACT)*( FEBMIN( IHIGH) +FEBMIN(IHIGH))/2.
             END IF
      IF(IPAR.EQ. 14) THEN
             FLCOSL= SOLACT*( FMMAX( ILOW ) +FMMAX(ILOW ))/2. +
     &          (1.-SOLACT)*( FMMIN( ILOW ) +FMMIN(ILOW ))/2.
             FLCOSH= SOLACT*( FMMAX( IHIGH) +FMMAX(IHIGH))/2. +
     &          (1.-SOLACT)*( FMMIN( IHIGH) +FMMIN(IHIGH))/2.
             END IF
      IF(IPAR.EQ.-14) THEN
             FLCOSL= SOLACT*( FMBMAX( ILOW ) +FMBMAX(ILOW ))/2. +
     &          (1.-SOLACT)*( FMBMIN( ILOW ) +FMBMIN(ILOW ))/2.
             FLCOSH= SOLACT*( FMBMAX( IHIGH) +FMBMAX(IHIGH))/2. +
     &          (1.-SOLACT)*( FMBMIN( IHIGH) +FMBMIN(IHIGH))/2.
             END IF
      IF(ILOW.NE.IHIGH) THEN
C       PRINT *,'I,ENELOW L/H FLCOS',
C    &        ILOW, IHIGH, ENELOW(ILOW),ENELOW(IHIGH),
C    &        (ENELOW(ILOW)/ENELOW(IHIGH)),FLCOSL,FLCOSH
        B     = ALOG (FLCOSL/FLCOSH) / ALOG (ENELOW(ILOW)/ENELOW(IHIGH))
        A     = FLCOSL / ENELOW(ILOW)**B
        FLUXNE= A*ENEU**B
      ELSE
        FLUXNE=FLCOSL
      END IF
C
C CONSIDER ANGULAR DISTRIBUTION USING THE CALCULATION FROM HONDA.
C 
      DDIR=SQRT(DIR(1)**2+DIR(2)**2+DIR(3)**2)
      IF(DDIR.EQ.0.) RETURN
*      IF(DDIR.NE.1.) THEN
*      DIR(1)=DIR(1)/DDIR
*      DIR(2)=DIR(2)/DDIR
*      DIR(3)=DIR(3)/DDIR
*      ENDIF
C
      FLUXHT=0.
      DO 3500 IDIR=1,20
      DIRCAL(3)=-1.05+0.1*FLOAT(IDIR)
      DIRCAL(1) = SQRT (1.-DIRCAL(3)**2)
      DIRCAL(2) = 0.
      CALL NEHONFLX(IPAR,ENEU,DIRCAL,FLUXHO)
      FLUXHT=FLUXHT+FLUXHO
 3500 CONTINUE
      FLUXHO=FLUXHT/20.
      RATIO=FLUXNE/FLUXHO
      CALL NEHONFLX(IPAR,ENEU,DIR,FLUXHO)
      FLUXNE=FLUXHO*RATIO
      GO TO 1000
C
C
C
C -- E(NEU) BETWEEN 3 GEV AND 10 GEV
  200 CONTINUE
*      DDIR=SQRT(DIR(1)**2+DIR(2)**2+DIR(3)**2)
*      IF(DDIR.NE.1. .AND. DDIR.NE.0.) THEN
*      DIR(1)=DIR(1)/DDIR
*      DIR(2)=DIR(2)/DDIR
*      DIR(3)=DIR(3)/DDIR
*      ENDIF
      COST = DIR(3)
      DO 220 ICOS=2,8
*      IF (ABS(COST).GT.PCOSH(ICOS)) THEN
      IF (ABS(COST).GE.PCOSH(ICOS)) THEN
                         ICOS3 = ICOS-1
                         ICOS4 = ICOS
                         DCOS3 = - ( ABS(COST)-PCOSH(ICOS-1) )
                         DCOS4 =     ABS(COST)-PCOSH(ICOS)
                         GO TO 221
                         END IF
 220  CONTINUE
      RETURN
 221  CONTINUE
C
      IF(IPAR.EQ. 12) THEN
             FLCOSL= SOLACT*( FEMAX( ILOW ) +FEMAX(ILOW ))/2. +
     &          (1.-SOLACT)*( FEMIN( ILOW ) +FEMIN(ILOW ))/2.
             FLCOSH= ( DCOS4*FEHI(ICOS3,IHIGH)+DCOS3*FEHI(ICOS4,IHIGH))
     &                      /(DCOS3+DCOS4)*FNENEB(COST)
             END IF
      IF(IPAR.EQ.-12) THEN
             FLCOSL= SOLACT*( FEBMAX( ILOW ) +FEBMAX(ILOW ))/2. +
     &          (1.-SOLACT)*( FEBMIN( ILOW ) +FEBMIN(ILOW ))/2.
             FLCOSH= ( DCOS4*FEHI(ICOS3,IHIGH)+DCOS3*FEHI(ICOS4,IHIGH))
     &                      /(DCOS3+DCOS4)*(1.-FNENEB(COST))
             END IF
      IF(IPAR.EQ. 14) THEN
             FLCOSL= SOLACT*( FMMAX( ILOW ) +FMMAX(ILOW ))/2. +
     &          (1.-SOLACT)*( FMMIN( ILOW ) +FMMIN(ILOW ))/2.
             FLCOSH= ( DCOS4*FMHI(ICOS3,IHIGH)+DCOS3*FMHI(ICOS4,IHIGH))
     &                      /(DCOS3+DCOS4)*FNMNMB(COST)
             END IF
      IF(IPAR.EQ.-14) THEN
             FLCOSL= SOLACT*( FMBMAX( ILOW ) +FMBMAX(ILOW ))/2. +
     &          (1.-SOLACT)*( FMBMIN( ILOW ) +FMBMIN(ILOW ))/2.
             FLCOSH= ( DCOS4*FMHI(ICOS3,IHIGH)+DCOS3*FMHI(ICOS4,IHIGH))
     &                      /(DCOS3+DCOS4)*(1.-FNMNMB(COST))
             END IF
C
C CONSIDER ANGULAR DISTRIBUTION USING THE CALCULATION FROM HONDA.
      FLUXHT=0.
      DO 3510 IDIR=1,20
      DIRCAL(3)=-1.05+0.1*FLOAT(IDIR)
      DIRCAL(1) = SQRT (1.-DIRCAL(3)**2)
      DIRCAL(2) = 0.
      CALL NEHONFLX ( IPAR,ENELOW(ILOW),DIRCAL,FLUXHO)
      FLUXHT=FLUXHT+FLUXHO
 3510 CONTINUE
      FLUXHO=FLUXHT/20.
      RATIO=FLCOSL/FLUXHO
      CALL NEHONFLX ( IPAR,ENELOW(ILOW),DIR,FLUXHO)
      FLCOSL=FLUXHO*RATIO
C
      B     = ALOG (FLCOSL/FLCOSH) / ALOG (ENELOW(ILOW)/ENEHI(IHIGH))
      A     = FLCOSL / ENELOW(ILOW)**B
      FLUXNE= A*ENEU**B
      GO TO 1000
C
C
C
C -- E(NEU) LARGER THAN 10 GEV
  300 CONTINUE
*      DDIR=SQRT(DIR(1)**2+DIR(2)**2+DIR(3)**2)
*      IF(DDIR.NE.1. .AND. DDIR.NE.0.) THEN
*      DIR(3)=DIR(3)/DDIR
*      ENDIF
      COST = DIR(3)
      DO 320 ICOS=2,8
*      IF (ABS(COST).GT.PCOSH(ICOS)) THEN
      IF (ABS(COST).GE.PCOSH(ICOS)) THEN
                         ICOS3 = ICOS-1
                         ICOS4 = ICOS
                         DCOS3 = - ( ABS(COST)-PCOSH(ICOS-1) )
                         DCOS4 =     ABS(COST)-PCOSH(ICOS)
                         GO TO 321
                         END IF
 320  CONTINUE
      RETURN
 321  CONTINUE
C
      IF(IPAR.EQ. 12) THEN
             FLCOSL= ( DCOS4*FEHI(ICOS3,ILOW )+DCOS3*FEHI(ICOS4,ILOW ))
     &                      /(DCOS3+DCOS4)*FNENEB(COST)
             FLCOSH= ( DCOS4*FEHI(ICOS3,IHIGH)+DCOS3*FEHI(ICOS4,IHIGH))
     &                      /(DCOS3+DCOS4)*FNENEB(COST)
             END IF
      IF(IPAR.EQ.-12) THEN
             FLCOSL= ( DCOS4*FEHI(ICOS3,ILOW )+DCOS3*FEHI(ICOS4,ILOW ))
     &                      /(DCOS3+DCOS4)*(1.-FNENEB(COST))
             FLCOSH= ( DCOS4*FEHI(ICOS3,IHIGH)+DCOS3*FEHI(ICOS4,IHIGH))
     &                      /(DCOS3+DCOS4)*(1.-FNENEB(COST))
             END IF
      IF(IPAR.EQ. 14) THEN
             FLCOSL= ( DCOS4*FMHI(ICOS3,ILOW )+DCOS3*FMHI(ICOS4,ILOW ))
     &                      /(DCOS3+DCOS4)*FNMNMB(COST)
             FLCOSH= ( DCOS4*FMHI(ICOS3,IHIGH)+DCOS3*FMHI(ICOS4,IHIGH))
     &                      /(DCOS3+DCOS4)*FNMNMB(COST)
             END IF
      IF(IPAR.EQ.-14) THEN
             FLCOSL= ( DCOS4*FMHI(ICOS3,ILOW )+DCOS3*FMHI(ICOS4,ILOW ))
     &                      /(DCOS3+DCOS4)*(1.-FNMNMB(COST))
             FLCOSH= ( DCOS4*FMHI(ICOS3,IHIGH)+DCOS3*FMHI(ICOS4,IHIGH))
     &                      /(DCOS3+DCOS4)*(1.-FNMNMB(COST))
             END IF
      IF(ILOW.NE.IHIGH) THEN
        B     = ALOG (FLCOSL/FLCOSH) / ALOG (ENEHI(ILOW)/ENEHI(IHIGH))
        A     = FLCOSL / ENEHI(ILOW)**B
        FLUXNE= A*ENEU**B
      ELSE
        FLUXNE=FLCOSL
      END IF
C
C -- SET NEUTRINO FLUX
C
 1000 CONTINUE
      FLUXNE=FLUXNE*1.0
C
      RETURN
      END

*******************************************************************
*     ------------------------------
      SUBROUTINE NEFLUX(IPAR,E,FLUX)
*     ------------------------------
*
*     ( Purpose )
*       GIVE NEUTRINO FLUX IN UNIT (EVENTS/M**2/SEC/STR/GEV)
* 
*     ( Input )
*        IPAR : NEUTRINO TYPE
*               12 : NUE
*              -12 : NUEBAR
*               14 : NUMU
*              -14 : NUMUBAR
*        E    : ENERGY OF NEUTRINO (UNIT : GEV)
*
*     ( Output )
*        FLUX : 
*
*     ( Creation Date and Auther )
*        1983.11.13 ; FIRST VERSION BY M.NAKAHATA
*
*        1985.01.09 ; M.NAKAHATA
*                     LESS THAN 200 MEV DATA WAS ADDED. 
*                     (CALCULATION BY T.K.GAISSER)
*
*        1987.10.07 ; M.TAKITA  7-OCT-1987
*                     THE AVERAGE FLUX OF THE SOLAR MAX AND MIN 
*                     IS INSTALLED AND FLUX IS FIT AS F=A*E**B 
*                     INSTEAD OF F=A*10**BE]
*
*        1988.11.22 ; T.KAJITA  
*                     ### DATA ARE UPDATED ####
*                     LATEST GAISSER-FLUX IS INSTALLED   (NOV-88)
*                     VOLKOVA FLUX IS INSTALLED FOR E(NEU)>10GEV
*                     SOLAR MIN,MAX OR ANY COMBINATION OF MIN/MAX IS AVAILABLE
*                     IF(SOLACT=1.) -- SOLAR MAX.
*                     IF(SOLACT=0.) -- SOLAR MIN.
*                     (  IF(SOLACT=0.5) -- AVERAGE OF SOLAR MIN/MAX )
*                     NEU/(NEU-BAR) RATIO IS INSTALLED ACCORDING 
*                     TO THE GAISSER'S TABLE
*                     ALSO NEU/(NEU-BAR) AT HIGH ENERGY IS ESTIMATED 
*                     FROM THE SAME TABLE
*                     CORRECTION OF NEU(E)/NEU(MU) RATIO AT LOW E (<100MEV) 
*                     IS ESTIMATED FROM THE HIGHER ENERGY CALCULATION
*     
*        2004.02.04 ; C.Saji
*                     ADD TAU MODE
*
*        2007.08.27 ; G.Mitsuka, HKKM06 is supported
*
**********************************************************************
#include "necardap.h"
#include "nesolact.h"
      DIMENSION ENEU(19),
     &          FEMAX(19),FMMAX(19),FEBMAX(19),FMBMAX(19),
     &          FEMIN(19),FMMIN(19),FEBMIN(19),FMBMIN(19)
C
      DATA ENEU/  0.01,    0.025,    0.035,     0.05,     0.07,
     &            0.09,    0.3  ,    0.5  ,     0.7 ,     0.9 ,
     &            1.1 ,    1.3  ,    1.5  ,     1.7 ,     1.9 ,
     &            3.0 ,    10.  ,     30. ,    100.   /
C
      DATA FEMAX/  5.01E3,  7.84E3,  7.55E3,  6.57E3,  4.92E3,
     &             3.97E3,  9.35E2,  3.76E2,  1.91E2, 10.95E1 ,
     &             6.90E1,  4.56E1,  3.10E1, 2.265E1, 1.675E1 ,
     &            5.206E0,1.022E-1,2.422E-3,3.699E-5  /
      DATA FEBMAX/ 5.01E3,  7.34E3,  7.15E3,  6.06E3,  4.77E3 ,
     &             3.77E3,  8.70E2, 3.415E2, 1.695E2,  9.60E1 ,
     &             5.80E1, 3.805E1,  2.61E1,  1.90E1,  1.42E1 ,
     &            4.137E0,0.812E-1,1.923E-3,2.936E-5  /
      DATA FMMAX/  7.5E3 ,  1.3E4 ,  1.3E4 , 1.175E4,  9.75E3 ,
     &             8.0E3 , 1.880E3,  7.50E2, 3.780E2, 2.155E2 ,
     &            1.345E2,  9.05E1,  6.30E1, 4.525E1, 3.355E1 ,
     &            1.105E1,3.142E-1,1.150E-2,3.135E-4  /
      DATA FMBMAX/ 7.5E3 ,  1.3E4 ,  1.3E4 , 1.175E4,  9.75E3 ,
     &             8.0E3 ,  1.88E3,  7.45E2, 3.725E2,  2.13E2 ,
     &            1.330E2,  8.90E1,  6.15E1, 4.405E1, 3.295E1 ,
     &            1.069E1,2.858E-1,1.046E-2,2.850E-4  /
C
      DATA FEMIN/  6.01E3,  9.57E3,  9.25E3,  8.01E3,  5.95E3 ,
     &             4.76E3, 10.25E2, 3.995E2, 1.995E2,  1.13E2 ,
     &             7.05E1, 4.655E1, 3.155E1, 2.295E1, 1.695E1 ,
     &            5.175E0,1.022E-1,2.422E-3,3.699E-5  /
      DATA FEBMIN/ 6.01E3,  8.95E3,  8.77E3,  7.39E3,  5.77E3 ,
     &             4.52E3,  9.40E2, 3.585E2, 1.755E2,  9.85E1 ,
     &             5.95E1,  3.86E1,  2.64E1,  1.92E1,  1.43E1 ,
     &            4.125E0,0.812E-1,1.923E-3,2.936E-5  /
      DATA FMMIN/  9.0E3 ,  1.59E4,  1.59E4, 1.433E4, 11.80E3 ,
     &             9.6E3 , 2.045E3,  7.90E2,  3.91E2,  2.21E2 ,
     &             1.37E2,  9.20E1,  6.35E1,  4.56E1,  3.37E1 ,
     &            1.105E1,3.142E-1,1.150E-2,3.135E-4  /
      DATA FMBMIN/ 9.0E3 ,  1.59E4,  1.59E4, 1.433E4, 11.80E3 ,
     &             9.6E3 , 2.045E3,  7.85E2,  3.86E2, 2.185E2 ,
     &             1.36E2,  9.00E1,  6.20E1,  4.44E1,  3.32E1 ,
     &            1.069E1,2.858E-1,1.046E-2,2.850E-4  /
C
      FLUX=0.
      IF(IABS(IPAR).NE.12 .AND. IABS(IPAR).NE.14
     $                    .AND. IABS(IPAR).NE.16 ) RETURN
      IF(IPAR.EQ.16) then
         IPAR_FLUX= 14
      ELSEIF(IPAR.EQ.-16) then
         IPAR_FLUX=-14
      ELSE
         IPAR_FLUX=IPAR
      ENDIF
C
      IF(ISOL.EQ.0) SOLACT=0.5
      IF(SOLACT.LT.0.) SOLACT=0.
      IF(SOLACT.GT.1.) SOLACT=1.

C--- CHOOSE FLUX
C---  NEFLXFLG = 0    : Honda06
C---             1    : Gaisser(96)
C---             2    : Lee    (Not supported now)
C---             3    : Naumov (Not supported now)
C---             4    : Gaisser(OLD)      
C---             5    : Honda (1D)      
C---             5    : Honda (1D)      
C---             6    : Honda03

C---HONDA---
      IF (NEFLXFLG.eq.0.or.NEFLXFLG.eq.6) THEN
         CALL NEFLXHON(IPAR_FLUX,E,FLUX)
         RETURN
      ENDIF
C---GAISSER(96)---
      IF (NEFLXFLG.eq.1) THEN
         CALL NEFLXGAI(IPAR_FLUX,E,FLUX)
         RETURN
      ENDIF
C---HONDA(1D)---
      IF (NEFLXFLG.eq.5) THEN
         CALL NEFLXHON_1D(IPAR_FLUX,E,FLUX)
         RETURN
      ENDIF
C---ERROR---
      IF (NEFLXFLG.eq.2.or.NEFLXFLG.eq.3.or.
     &     NEFLXFLG.eq.4) THEN
         write(*,*) "NEFLXFLG = ",NEFLXFLG," is not supported now."
         STOP
      ENDIF
C---GAISSER(OLD)---
C
C Change Y.H / June 3
C      IF(E.LE.0.01) RETURN
      IF(E.GE.100.) RETURN
 
      IF(E.LT.0.01) RETURN
C      IF(E.GT.100.) RETURN

C
      I=1
  100 IF(ENEU(I).GT.E)GO TO 200
      I=I+1
      GO TO 100
C
  200 CONTINUE
C
C--FIT FLUX BY F=A*E**B
C+NEU(E)
      IF(IPAR_FLUX.EQ. 12) THEN
         BMIN=ALOG(FEMIN(I-1)/FEMIN(I))/ALOG(ENEU(I-1)/ENEU(I))
         BMAX=ALOG(FEMAX(I-1)/FEMAX(I))/ALOG(ENEU(I-1)/ENEU(I))
         AMIN=FEMIN(I-1)/ENEU(I-1)**BMIN
         AMAX=FEMAX(I-1)/ENEU(I-1)**BMAX
      END IF
C+NEU(E)-BAR
      IF(IPAR_FLUX.EQ.-12) THEN
         BMIN=ALOG(FEBMIN(I-1)/FEBMIN(I))/ALOG(ENEU(I-1)/ENEU(I))
         BMAX=ALOG(FEBMAX(I-1)/FEBMAX(I))/ALOG(ENEU(I-1)/ENEU(I))
         AMIN=FEBMIN(I-1)/ENEU(I-1)**BMIN
         AMAX=FEBMAX(I-1)/ENEU(I-1)**BMAX
      END IF
C+NEU(MU)
      IF(IPAR_FLUX.EQ. 14) THEN
         BMIN=ALOG(FMMIN(I-1)/FMMIN(I))/ALOG(ENEU(I-1)/ENEU(I))
         BMAX=ALOG(FMMAX(I-1)/FMMAX(I))/ALOG(ENEU(I-1)/ENEU(I))
         AMIN=FMMIN(I-1)/ENEU(I-1)**BMIN
         AMAX=FMMAX(I-1)/ENEU(I-1)**BMAX
      END IF
C+NEU(MU)-BAR
      IF(IPAR_FLUX.EQ.-14) THEN
         BMIN=ALOG(FMBMIN(I-1)/FMBMIN(I))/ALOG(ENEU(I-1)/ENEU(I))
         BMAX=ALOG(FMBMAX(I-1)/FMBMAX(I))/ALOG(ENEU(I-1)/ENEU(I))
         AMIN=FMBMIN(I-1)/ENEU(I-1)**BMIN
         AMAX=FMBMAX(I-1)/ENEU(I-1)**BMAX
      END IF
C
      FLUX=SOLACT*(AMAX*E**BMAX)+(1.-SOLACT)*(AMIN*E**BMIN)
      IF(E.GT.2.) RETURN
C
C
C +++ SMALL CORRECTION TO THE CALCULATED FLUX
C    THIS IS BECAUSE THE GAISSERS TABLE IS THE INTEGRATED FLUX
C    BETWEEN E AND E+DE WITH DE=200MEV, AND BECAUSE WE SET THE NEUTRINO
C    ENERGY TO BE E+DE/2. THEREFORE THE ABOVE FLUX OVERESTIMATES BY
C    1-2%.
C
      FLUX=FLUX*(1.-0.02*(2.-E))
C
      RETURN
      END

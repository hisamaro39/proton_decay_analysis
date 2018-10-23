***********************************************************************
*     ----------------------------------
      SUBROUTINE NENAUFLX(IPAR,E0,FLUXN)
*     ----------------------------------
*
*     ( Purpose )
*       FUNCTION TO GIVE NEUTRINO FLUX BASED ON NAUMOV FLUX
*
*     ( Input )
*       IPAR : NEUTRINO TYPE
*             12 : NUE
*            -12 : NUEBAR
*             14 : NUMU
*            -14 : NUMUBAR
*       E0   : ENERGY OF NEUTRINO (UNIT : GEV)
*     
*     ( Output )
*       FLUXN: FLUX OF NEUTRINO FOR GIVEN IPAR AND E 
*                                   (UNIT : EVENTS/M**2/SEC/STR/GEV)
*
*     ( Creation Date and Auther )
*       1992.09.03 ; T.KAJITA  
*                    VOLKOVA FLUX IS INSTALLED FOR E(NEU)>10GEV
*       1996.01.23 ; Y.KANAYA
*                    FOR SK   
*
*     ( Comment )
*       SOLAR MIN,MAX OR ANY COMBINATION OF MIN/MAX IS NOT AVAILABLE
*            IF(SOLACT=1.) -- SOLAR MAX.
*            IF(SOLACT=0.) -- SOLAR MIN.
*         (  IF(SOLACT=0.5) -- AVERAGE OF SOLAR MIN/MAX )
*       CALL NERDNAU BEFOREHAND
*
*********************************************************************
#include "nesolact.h"
      COMMON /NENAUMF / ENEU(33),
     &       FEMIN(33), FEMAX(33), FEBMIN(33), FEBMAX(33),
     &       FMMIN(33), FMMAX(33), FMBMIN(33), FMBMAX(33)
      DATA ICALL/0/
C
C --  READ TABLE
C
      IF (ICALL.EQ.0) THEN
         ICALL = 1
         CALL NERDNAU
      ENDIF
      
C     
      IF(IABS(IPAR).NE.12 .AND. IABS(IPAR).NE.14)RETURN
      IF(E0.LE.0.01) RETURN
      IF(E0.GE.100.) RETURN
C
      IF(ISOL.EQ.0) SOLACT=0.5
      IF(SOLACT.LT.0.) SOLACT=0.
      IF(SOLACT.GT.1.) SOLACT=1.
C
      FLUXN=0.
C
      IF ( E0.GT.10. ) GO TO 1000
C
      E=E0
      IF ( E.GT.ENEU(33) ) E=ENEU(33)
      IF ( E.LE.ENEU(1) ) E=ENEU(1)+0.0001
C
      I=2
  100 IF(ENEU(I).GT.E)GO TO 200
      I=I+1
        IF(I.GE.34) THEN
        I=33
        GO TO 200
        END IF
      GO TO 100
C
C--FIT FLUX BY F=A*E**B
  200 CONTINUE
C+NEU(E)
      IF(IPAR.EQ. 12) THEN
         BMIN=ALOG(FEMIN(I-1)/FEMIN(I))/ALOG(ENEU(I-1)/ENEU(I))
         BMAX=ALOG(FEMAX(I-1)/FEMAX(I))/ALOG(ENEU(I-1)/ENEU(I))
         AMIN=FEMIN(I-1)/ENEU(I-1)**BMIN
         AMAX=FEMAX(I-1)/ENEU(I-1)**BMAX
      END IF
C+NEU(E)-BAR
      IF(IPAR.EQ.-12) THEN
         BMIN=ALOG(FEBMIN(I-1)/FEBMIN(I))/ALOG(ENEU(I-1)/ENEU(I))
         BMAX=ALOG(FEBMAX(I-1)/FEBMAX(I))/ALOG(ENEU(I-1)/ENEU(I))
         AMIN=FEBMIN(I-1)/ENEU(I-1)**BMIN
         AMAX=FEBMAX(I-1)/ENEU(I-1)**BMAX
      END IF
C+NEU(MU)
      IF(IPAR.EQ. 14) THEN
         BMIN=ALOG(FMMIN(I-1)/FMMIN(I))/ALOG(ENEU(I-1)/ENEU(I))
         BMAX=ALOG(FMMAX(I-1)/FMMAX(I))/ALOG(ENEU(I-1)/ENEU(I))
         AMIN=FMMIN(I-1)/ENEU(I-1)**BMIN
         AMAX=FMMAX(I-1)/ENEU(I-1)**BMAX
      END IF
C+NEU(MU)-BAR
      IF(IPAR.EQ.-14) THEN
         BMIN=ALOG(FMBMIN(I-1)/FMBMIN(I))/ALOG(ENEU(I-1)/ENEU(I))
         BMAX=ALOG(FMBMAX(I-1)/FMBMAX(I))/ALOG(ENEU(I-1)/ENEU(I))
         AMIN=FMBMIN(I-1)/ENEU(I-1)**BMIN
         AMAX=FMBMAX(I-1)/ENEU(I-1)**BMAX
      END IF
C
      FLUXN=SOLACT*(AMAX*E**BMAX)+(1.-SOLACT)*(AMIN*E**BMIN)
      IF ( E0.GT.ENEU(33) ) GOTO 1000
      RETURN
C
 1000 CONTINUE
      FLUXN3=FLUXN
      CALL NEFLUX(IPAR,E0,FLUXN)
      IF(E0.GT.10.) RETURN
      CALL NEFLUX(IPAR,ENEU(33),FLUXG3)
      FLUXN = FLUXN *(1.- (10.-E0)/(10.-ENEU(33))*(1.-FLUXN3/FLUXG3) )
C
      RETURN
      END


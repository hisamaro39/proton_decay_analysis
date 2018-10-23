***********************************************************************
*     -----------------------------------------
      SUBROUTINE NEHONFLX(IPAR,ENEU,DIR,FLUXNE)
*     -----------------------------------------
*
*     ( Purpose )
*       GET THE ATMOSPHERIC NEUTRINO FLUX FOR A GIVEN NENU,DIR,IPAR
*       THE FLUX IS THE ONE CALCULATED BY HONDA IN 1990 (LOW E) AND
*       1993 (HIGH E).
*
*     ( Input )
*       IPAR; NEUTRINO TYPE
*             12 ; NUE
*            -12 ; NUEBAR
*             14 ; NUMU
*            -14 ; NUMUBAR
*       ENEU; NEUTRINO ENERGY IN GEV
*       DIR ; NEUTRINO DIRECTION
*
*     ( Output )
*       FLUXNE; NEUTRINO FLUX   ( /M**2/SEC/SR/GEV )
*
*     ( Creation Date and Auther )
*       1993.10.22 ; UPDATED BY T.KAJTIA
*                    LATEST LOW-E  HONDA-FLUX IS INSTALLED (<4GEV) (1990)
*                    LATEST HIGH-E HONDA-FLUX IS INSTALLED (>5GEV) (1993)
*                    NEU/(NEU-BAR) RATIO IS INSTALLED ACCORDING TO THE 
*                    HONDA'S TABLE
*       1996.01.24 ; Y.KANAYA
*                    FOR SK
*
*     ( Comment )
*       SOLAR MIN,MAX OR ANY COMBINATION OF MIN/MAX IS NOT AVAILABLE
*            IF(SOLACT=1.) -- SOLAR MAX.                  =MEAN
*            IF(SOLACT=0.) -- SOLAR MIN.                  =MEAN
*         (  IF(SOLACT=0.5) -- AVERAGE OF SOLAR MIN/MAX ) =MEAN
*       CALL NERDHON BEFOREHAND
*
*************************************************************************
#include "nesolact.h"
#include "skcnst.h"
C           NUMBER OF LOW- AND HIGH-ENERGY FLUX'S ENERGY BINS
      PARAMETER (NL=73,NH=41)
C
      DIMENSION DIR(3)
      COMMON /NEHONDAF/ ENELOW(NL), PCOSL(20),
     &       FEMIN(20,NL), FEMAX(20,NL), FEBMIN(20,NL), FEBMAX(20,NL),
     &       FMMIN(20,NL), FMMAX(20,NL), FMBMIN(20,NL), FMBMAX(20,NL)
      COMMON /NEHONDAH/ ENEHI(NH), PCOSH(20),
     &       FEHMN(20,NH), FEHMX(20,NH), FEBHMN(20,NH), FEBHMX(20,NH),
     &       FMHMN(20,NH), FMHMX(20,NH), FMBHMN(20,NH), FMBHMX(20,NH)
C
C         NEUTRINO ENERGY CUT FOR USING LOW-E FILE (ECUTL) AND
C         HIGH-E FILE (ECUTH).      UNIT=GEV.
      DATA ECUTL/3.5/, ECUTH/5.0/
      DATA ICALL/0/
C
C
C
      WRITE(0,*)"************************************************"
      WRITE(0,*)"*                                              *" 
      WRITE(0,*)"*                                              *"
      WRITE(0,*)"*  THIS SUBROUTINE(NEHONFLX) IS OBSOLETE!!!    *"
      WRITE(0,*)"*                                              *"
      WRITE(0,*)"*    EXECUTION TERMINATED                      *"
      WRITE(0,*)"*                                              *"
      WRITE(0,*)"*                                              *" 
      WRITE(0,*)"************************************************"
      STOP

C
C --  READ TABLE
C
      IF (ICALL.EQ.0) THEN
         ICALL = 1
         CALL NERDHON
      ENDIF

C
C -- CHECK NERDHON
      FLUXNL=0.
      FLUXNH=0.
      FLUXNE=0.
      IENERG=0
C --SOLAR ACTIVITY
      IF(ISOL.EQ.0) SOLACT=0.5
      IF(SOLACT.LT.0.) SOLACT=0.
      IF(SOLACT.GT.1.) SOLACT=1.
C
C
C --DETERMINE ZENITH-ANGLE BIN
      DDIR=SQRT(DIR(1)**2+DIR(2)**2+DIR(3)**2)
      IF(DDIR.LE.0.) THEN
      CALL NEFINT(ENEU,IPAR,4,FLUXNE)
      FLUXNE=FLUXNE/4/PI
      RETURN
*      ELSEIF(DDIR.NE.1.) THEN
*      DIR(1)=DIR(1)/DDIR
*      DIR(2)=DIR(2)/DDIR
*      DIR(3)=DIR(3)/DDIR
      ENDIF
C
      COST = DIR(3)
      DO 110 ICOS=2,20
       IF (COST.GT.PCOSL(ICOS)) THEN
*       IF (COST.GE.PCOSL(ICOS)) THEN
           ICOS1 = ICOS-1
           ICOS2 = ICOS
           DCOS1 = - ( COST-PCOSL(ICOS-1) )
           DCOS2 =     COST-PCOSL(ICOS)
           GO TO 111
       END IF
 110  CONTINUE
      RETURN
 111  CONTINUE
C
C
C -- E(NEU) LESS THAN 7 GEV (LOW-ENERGY FILE)
C
C --DETERMINE ENERGY BIN
C
      DO 1 I=1,NL
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
      GO TO 1000
C
C -- CALCULATE FLUX
C
   10 CONTINUE
      IF(IPAR.EQ. 12) THEN
             FLCOSL= SOLACT*( DCOS2*FEMAX(ICOS1,ILOW )
     &                       +DCOS1*FEMAX(ICOS2,ILOW ))/(DCOS1+DCOS2) +
     &          (1.-SOLACT)*( DCOS2*FEMIN(ICOS1,ILOW )
     &                       +DCOS1*FEMIN(ICOS2,ILOW ))/(DCOS1+DCOS2)
             FLCOSH= SOLACT*( DCOS2*FEMAX(ICOS1,IHIGH)
     &                       +DCOS1*FEMAX(ICOS2,IHIGH))/(DCOS1+DCOS2) +
     &          (1.-SOLACT)*( DCOS2*FEMIN(ICOS1,IHIGH)
     &                       +DCOS1*FEMIN(ICOS2,IHIGH))/(DCOS1+DCOS2)
             END IF
      IF(IPAR.EQ.-12) THEN
             FLCOSL= SOLACT*( DCOS2*FEBMAX(ICOS1,ILOW )
     &                       +DCOS1*FEBMAX(ICOS2,ILOW ))/(DCOS1+DCOS2) +
     &          (1.-SOLACT)*( DCOS2*FEBMIN(ICOS1,ILOW )
     &                       +DCOS1*FEBMIN(ICOS2,ILOW ))/(DCOS1+DCOS2)
             FLCOSH= SOLACT*( DCOS2*FEBMAX(ICOS1,IHIGH)
     &                       +DCOS1*FEBMAX(ICOS2,IHIGH))/(DCOS1+DCOS2) +
     &          (1.-SOLACT)*( DCOS2*FEBMIN(ICOS1,IHIGH)
     &                       +DCOS1*FEBMIN(ICOS2,IHIGH))/(DCOS1+DCOS2)
             END IF
      IF(IPAR.EQ. 14) THEN
             FLCOSL= SOLACT*( DCOS2*FMMAX(ICOS1,ILOW )
     &                       +DCOS1*FMMAX(ICOS2,ILOW ))/(DCOS1+DCOS2) +
     &          (1.-SOLACT)*( DCOS2*FMMIN(ICOS1,ILOW )
     &                       +DCOS1*FMMIN(ICOS2,ILOW ))/(DCOS1+DCOS2)
             FLCOSH= SOLACT*( DCOS2*FMMAX(ICOS1,IHIGH)
     &                       +DCOS1*FMMAX(ICOS2,IHIGH))/(DCOS1+DCOS2) +
     &          (1.-SOLACT)*( DCOS2*FMMIN(ICOS1,IHIGH)
     &                       +DCOS1*FMMIN(ICOS2,IHIGH))/(DCOS1+DCOS2)
             END IF
      IF(IPAR.EQ.-14) THEN
             FLCOSL= SOLACT*( DCOS2*FMBMAX(ICOS1,ILOW )
     &                       +DCOS1*FMBMAX(ICOS2,ILOW ))/(DCOS1+DCOS2) +
     &          (1.-SOLACT)*( DCOS2*FMBMIN(ICOS1,ILOW )
     &                       +DCOS1*FMBMIN(ICOS2,ILOW ))/(DCOS1+DCOS2)
             FLCOSH= SOLACT*( DCOS2*FMBMAX(ICOS1,IHIGH)
     &                       +DCOS1*FMBMAX(ICOS2,IHIGH))/(DCOS1+DCOS2) +
     &          (1.-SOLACT)*( DCOS2*FMBMIN(ICOS1,IHIGH)
     &                       +DCOS1*FMBMIN(ICOS2,IHIGH))/(DCOS1+DCOS2)
             END IF
C
      IF(ILOW.NE.IHIGH) THEN
        B     = ALOG (FLCOSL/FLCOSH) / ALOG (ENELOW(ILOW)/ENELOW(IHIGH))
        A     = FLCOSL / ENELOW(ILOW)**B
        FLUXNL= A*ENEU**B
      ELSE
        FLUXNL=FLCOSL
      END IF
C
C
C
C -- E(NEU) LARGER THAN 1 GEV (HIGH-ENERGY FILE)
C
C --DETERMINE ENERGY BIN
C
 1000 CONTINUE
      IF(ENEU.LT.ENEHI(1)) GO TO 2000
      DO 1001 I=2,NH
      IMAX=I
      DENEU = ENEU-ENEHI(I)
       IF(DENEU.LT.0.) THEN
           ILOW=I-1
           IHIGH=I
           IENERG=1
           GO TO 210
       END IF
 1001 CONTINUE
      RETURN
C
C -- CALCULATE FLUX
C
  210 CONTINUE
      IF(IPAR.EQ. 12) THEN
             FLCOSL= SOLACT*( DCOS2*FEHMX(ICOS1,ILOW )
     &                       +DCOS1*FEHMX(ICOS2,ILOW ))/(DCOS1+DCOS2) +
     &          (1.-SOLACT)*( DCOS2*FEHMN(ICOS1,ILOW )
     &                       +DCOS1*FEHMN(ICOS2,ILOW ))/(DCOS1+DCOS2)
             FLCOSH= SOLACT*( DCOS2*FEHMX(ICOS1,IHIGH)
     &                       +DCOS1*FEHMX(ICOS2,IHIGH))/(DCOS1+DCOS2) +
     &          (1.-SOLACT)*( DCOS2*FEHMN(ICOS1,IHIGH)
     &                       +DCOS1*FEHMN(ICOS2,IHIGH))/(DCOS1+DCOS2)
             END IF
      IF(IPAR.EQ.-12) THEN
             FLCOSL= SOLACT*( DCOS2*FEBHMX(ICOS1,ILOW )
     &                       +DCOS1*FEBHMX(ICOS2,ILOW ))/(DCOS1+DCOS2) +
     &          (1.-SOLACT)*( DCOS2*FEBHMN(ICOS1,ILOW )
     &                       +DCOS1*FEBHMN(ICOS2,ILOW ))/(DCOS1+DCOS2)
             FLCOSH= SOLACT*( DCOS2*FEBHMX(ICOS1,IHIGH)
     &                       +DCOS1*FEBHMX(ICOS2,IHIGH))/(DCOS1+DCOS2) +
     &          (1.-SOLACT)*( DCOS2*FEBHMN(ICOS1,IHIGH)
     &                       +DCOS1*FEBHMN(ICOS2,IHIGH))/(DCOS1+DCOS2)
             END IF
      IF(IPAR.EQ. 14) THEN
             FLCOSL= SOLACT*( DCOS2*FMHMX(ICOS1,ILOW )
     &                       +DCOS1*FMHMX(ICOS2,ILOW ))/(DCOS1+DCOS2) +
     &          (1.-SOLACT)*( DCOS2*FMHMN(ICOS1,ILOW )
     &                       +DCOS1*FMHMN(ICOS2,ILOW ))/(DCOS1+DCOS2)
             FLCOSH= SOLACT*( DCOS2*FMHMX(ICOS1,IHIGH)
     &                       +DCOS1*FMHMX(ICOS2,IHIGH))/(DCOS1+DCOS2) +
     &          (1.-SOLACT)*( DCOS2*FMHMN(ICOS1,IHIGH)
     &                       +DCOS1*FMHMN(ICOS2,IHIGH))/(DCOS1+DCOS2)
             END IF
      IF(IPAR.EQ.-14) THEN
             FLCOSL= SOLACT*( DCOS2*FMBHMX(ICOS1,ILOW )
     &                       +DCOS1*FMBHMX(ICOS2,ILOW ))/(DCOS1+DCOS2) +
     &          (1.-SOLACT)*( DCOS2*FMBHMN(ICOS1,ILOW )
     &                       +DCOS1*FMBHMN(ICOS2,ILOW ))/(DCOS1+DCOS2)
             FLCOSH= SOLACT*( DCOS2*FMBHMX(ICOS1,IHIGH)
     &                       +DCOS1*FMBHMX(ICOS2,IHIGH))/(DCOS1+DCOS2) +
     &          (1.-SOLACT)*( DCOS2*FMBHMN(ICOS1,IHIGH)
     &                       +DCOS1*FMBHMN(ICOS2,IHIGH))/(DCOS1+DCOS2)
             END IF
C
        B     = ALOG (FLCOSL/FLCOSH) / ALOG (ENEHI(ILOW)/ENEHI(IHIGH))
        A     = FLCOSL / ENEHI(ILOW)**B
        FLUXNH= A*ENEU**B
C
C
C
C -- SET NEUTRINO FLUX
C
 2000 CONTINUE
      IF(ENEU.LT.ECUTL) THEN
                   FLUXNE=FLUXNL
      ELSE IF(ENEU.LT.ECUTH) THEN
                   DECUT=ECUTH-ECUTL
                   DENEU=ENEU-ECUTL
                   FLUXNE=((DECUT-DENEU)*FLUXNL+DENEU*FLUXNH)/DECUT
      ELSE
                   FLUXNE=FLUXNH
      END IF
C
      FLUXNE=FLUXNE*1.0
C
      RETURN
      END

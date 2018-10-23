***********************************************************************
*     --------------------------------------
      SUBROUTINE NEFINT(IPAR,ENEU,IFLG,FINT)
*     --------------------------------------
*
*     ( Purpose )
*       GET THE ATMOSPHERIC NEUTRINO FLUX FOR A GIVEN NENU,IPAR
*
*     ( Input )
*       IPAR; NEUTRINO TYPE
*             12 ; NUE
*            -12 ; NUEBAR
*             14 ; NUMU
*            -14 ; NUMUBAR
*       ENEU; NEUTRINO ENERGY IN GEV
*       IFLG; 
*              1   ; NEGAIFLX
*              2   ; NENAUFLX
*              3   ; NELEEFLX
*              ELSE; NEHONFLX
*
*     ( Output )
*       FINT; INITEGRATED NEUTRINO FLUX (/M**2/SEC/GEV)
*
*     ( Creation Date and Auther )
*       1996.01.25 :
*
*     ( Comment )
*
*************************************************************************
#include "nesolact.h"
#include "skcnst.h"
C
      DIMENSION DIRCAL(3)
C
C --  NEGAIFLX
C
      IF (IFLG.EQ.1) THEN
         FLUXHT = 0.0
         DO 10 IDIR = 1,20
            DIRCAL(3)=-1.05+0.1*FLOAT(IDIR)
            DIRCAL(1) = SQRT (1.-DIRCAL(3)**2)
            DIRCAL(2) = 0.
            CALL NEGAIFLX(IPAR,ENEU,DIRCAL,FLUXHO)
            FLUXHT = FLUXHT + FLUXHO * 0.1
 10      CONTINUE
         FINT = FLUXHT * 2.0 * PI
C     
C -- NENAUFLX
C      
      ELSE IF (IFLG.EQ.2) THEN
         CALL NENAUFLX(IPAR,ENEU,FINT)
         FINT=FINT*4.*PI
C 
C -- NELEEFLX 
C      
      ELSE IF (IFLG.EQ.3) THEN
C
C --     E < 3 GeV
C         
         IF (ENEU.LT.3.) THEN
            DIRCAL(1)=0.
            DIRCAL(2)=0.
            DIRCAL(3)=0.
            CALL NELEEFLX(IPAR,ENEU,DIRCAL,FINT)
            FINT=FINT*4.*PI
         ELSE
C
C --     E > 3 GeV
C            
            FLUXHT = 0.0
            DO 20 IDIR = 1,20
               DIRCAL(3)=-1.05+0.1*FLOAT(IDIR)
               DIRCAL(1) = SQRT (1.-DIRCAL(3)**2)
               DIRCAL(2) = 0.
               CALL NELEEFLX(IPAR,ENEU,DIRCAL,FLUXHO)
               FLUXHT = FLUXHT + FLUXHO * 0.1
 20         CONTINUE
            FINT = FLUXHT * 2.0 * PI
         ENDIF
C     
C -- NEHONFLX
C            
      ELSE
         FLUXHT = 0.0
         DO 30 IDIR = 1,20
            DIRCAL(3)=-1.05+0.1*FLOAT(IDIR)
            DIRCAL(1) = SQRT (1.-DIRCAL(3)**2)
            DIRCAL(2) = 0.
            CALL NEHONFLX(IPAR,ENEU,DIRCAL,FLUXHO)
            FLUXHT = FLUXHT + FLUXHO * 0.1
 30      CONTINUE
         FINT = FLUXHT * 2.0 * PI
      ENDIF

      RETURN
      END

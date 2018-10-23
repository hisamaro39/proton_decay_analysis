************************************************************************
*     ---------------------------------
      SUBROUTINE neflxhon_1d(IPAR,ENEU,FLUX)
*     ---------------------------------
C     
*     (Purpose)
C        Honda Flux
C
C     (Input)
C       ENEU   : Energy of Neutrino(in GeV)
C       IPAR   : Flavor of neutrino
C
*     (Output)
C       FLUX   : Flux
C      
C     (Creation Date and Author)
C       1996.??.?? ; Y.Hayato
C       1997.??.?? ; Y.Hayato ( debug : Direction is changed to the
C                                       detector coordinates)       
C                                             
************************************************************************

      IMPLICIT NONE

#include "hondaflx_1d.h"
#include "nesolact.h"      
      
      INTEGER*4  IPAR
      REAL*4     ENEU,FLUX

      INTEGER*4  IPKIND

      INTEGER*4 ILOADED
      SAVE      ILOADED
      DATA      ILOADED/0/
      INTEGER*4 I
      REAL*4    BMIN,BMAX,AMIN,AMAX
      

C---- CHECK CALLED TWICE?
      IF (ILOADED.EQ.0) THEN
         CALL NERDHON_1D
         ILOADED = 1
      ENDIF

      IPKIND=0
      IF (IPAR.eq.12)  IPKIND=1
      IF (IPAR.eq.-12) IPKIND=2
      IF (IPAR.eq.14)  IPKIND=3
      IF (IPAR.eq.-14) IPKIND=4

      IF (IPKIND.eq.0) THEN
         write(*,*) "*** ERROR ***"
         write(*,*) "NEFLXHON:ILLEGAL IPAR=",IPAR
         STOP
      ENDIF

      FLUX=0.
      
      IF (ENEU.le.HONE_1D(1))   RETURN
      IF (ENEU.ge.HONE_1D(100)) RETURN

      I=1
 100  CONTINUE
         IF (HONE_1D(I).GT.ENEU) GOTO 200
         I=I+1
      GO TO 100

 200  CONTINUE

      BMIN=  ALOG(HFLX_1D(I-1,IPKIND,1)/HFLX_1D(I,IPKIND,1))
     $      /ALOG(HONE_1D(I-1)/HONE_1D(I))
      BMAX=  ALOG(HFLX_1D(I-1,IPKIND,2)/HFLX_1D(I,IPKIND,2))
     $      /ALOG(HONE_1D(I-1)/HONE_1D(I))
      AMIN= HFLX_1D(I-1,IPKIND,1)/HONE_1D(I-1)**BMIN
      AMAX= HFLX_1D(I-1,IPKIND,2)/HONE_1D(I-1)**BMAX

      FLUX=     SOLACT *(AMAX*ENEU**BMAX)
     $     +(1.-SOLACT)*(AMIN*ENEU**BMIN)

      RETURN
      END
      
      
      

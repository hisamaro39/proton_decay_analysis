************************************************************************
*     ---------------------------------
      SUBROUTINE neflxhon(IPAR,ENEU,FLUX)
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
C       2001.09.18 ; M.Ishitsuka 
C               replace flux table
C                  honda97mid.dat & honda96high.dat -> hkkm2.dat
C                                             
C       2007.12.05 ; G.Mitsuka, FC/PC's binning is merged with upmu
C
************************************************************************

      IMPLICIT NONE

#include "hondaflx.h"
#include "nesolact.h"      
      
      INTEGER*4  IPAR
      REAL*4     ENEU,FLUX

      INTEGER*4  IPKIND

      INTEGER*4 ILOADED
      SAVE      ILOADED
      DATA      ILOADED/0/
      INTEGER*4 I
      REAL*8    BMIN,BMAX,AMIN,AMAX

C---- CHECK CALLED TWICE?
      IF (ILOADED.EQ.0) THEN
         CALL NERDHON
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
      
      IF (ENEU.le.HONE(1))   RETURN
      IF (ENEU.ge.HONE(141)) RETURN

      I=1
 100  CONTINUE
         IF (HONE(I).GT.ENEU) GOTO 200
         I=I+1
      GO TO 100

 200  CONTINUE
      BMIN=  dble(ALOG(HFLX(I-1,IPKIND,1)/HFLX(I,IPKIND,1))
     $      /ALOG(HONE(I-1)/HONE(I)))
      BMAX=  dble(ALOG(HFLX(I-1,IPKIND,2)/HFLX(I,IPKIND,2))
     $      /ALOG(HONE(I-1)/HONE(I)))
      AMIN= dble(HFLX(I-1,IPKIND,1)/HONE(I-1)**BMIN)
      AMAX= dble(HFLX(I-1,IPKIND,2)/HONE(I-1)**BMAX)

      FLUX=     real(SOLACT *(AMAX*ENEU**BMAX)
     $     +(1.-SOLACT)*(AMIN*ENEU**BMIN))
      RETURN
      END

      
      
      

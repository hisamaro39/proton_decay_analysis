************************************************************************
*     ---------------------------------
      FUNCTION FNHONFX(ENEU,ORGDIR,IPAR)
*     ---------------------------------
C     
*     (Purpose)
C        Honda Differential Flux 
C
C     (Input)
C       ENEU   : Energy of Neutrino(in GeV)
C       ORGDIR : Direction of the neutrino
c             orgdir(3)<0 : downward direction
c             orgdir(3)>0 : upward direction
C       IPAR   : Flavor of neutrino
C
*     (Output)
C       Flux
C      
C     (Creation Date and Author)
C       1996.??.?? ; Y.Hayato
C       1997.??.?? ; Y.Hayato ( debug : Direction is changed to the
C                                       detector coordinates)       
C                                             
************************************************************************

C
      IMPLICIT NONE
#include"hondaflx.h"
#include"nesolact.h"
      
      REAL*4    FNHONFX
      REAL*4    ENEU
      REAL*4    ORGDIR(3)
      INTEGER*4 IPAR

      REAL*4    FHONL,FHONM

      REAL*4    FNHDFXL,FNHDFXM
      EXTERNAL  FNHDFXL,FNHDFXM

      REAL*4    DIR(3)

      INTEGER*4 ISTART
      DATA      ISTART/0/

      IF (ISTART.eq.0) THEN
         ISTART=1
         CALL NERDHON
      ENDIF

C---CHECK COMMON VAR(ISOL == 0 <--> SOLACT = 0.5)
C---AT THIS MOMENT THIS VALUE is MEANINGLESS..
      IF (ISOL.eq.0) SOLACT=0.50000

C--- GEOMETRY CONVERSION ---
C--- Invert Z -> -Z
C--- rotate 40.583 degree in X-Y plane

      DIR(3)=-1.*ORGDIR(3)
      DIR(1)= -(COS(3.141593*2435./60./180.)*ORGDIR(1)
     $       +SIN(3.141593*2435./60./180.)*ORGDIR(2))
      DIR(2)=-(-SIN(3.141593*2435./60./180.)*ORGDIR(1)
     $       +COS(3.141593*2435./60./180.)*ORGDIR(2))
    

      IF (ENEU.le.HONEM(1))  THEN
         FHONL=FNHDFXL(ENEU,DIR,IPAR)
         FNHONFX=FHONL
         RETURN
      ENDIF
      
      FHONM=FNHDFXM(ENEU,DIR,IPAR)
      FNHONFX=FHONM
      RETURN

      END

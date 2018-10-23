************************************************************************
*     ---------------------------------
      FUNCTION FNFLUFX(ENEU,ORGDIR,IPAR)
*     ---------------------------------
C     
*     (Purpose)
C        Fluka 3D Differential Flux 
C
C     (Input)
C       ENEU   : Energy of Neutrino(in GeV)
C       ORGDIR : Direction of the neutrino
C       IPAR   : Flavor of neutrino
C
*     (Output)
C       Flux
C      
C     (Creation Date and Author)
C       2000.10.14 ; A.Kibayashi
C                                             
************************************************************************

C
      IMPLICIT NONE
#include"flukaflx.h"
#include"nesolact.h"
      
      REAL*4    FNFLUFX
      REAL*4    ENEU
      REAL*4    ORGDIR(3)
      INTEGER*4 IPAR

      REAL*4    FHONL
      REAL*4    FFLU, FFLUH, FFLXHNRM

      REAL*4    fnhonfx
      EXTERNAL  fnhonfx

      REAL*4    FNFDFX
      EXTERNAL  FNFDFX

      REAL*4    DIR(3)

      INTEGER*4 ISTART
      DATA      ISTART/0/

      IF (ISTART.eq.0) THEN
         ISTART=1
         CALL NERDFLU
      ENDIF

C---CHECK COMMON VAR(ISOL == 0 <--> SOLACT = 0.5)
C---AT THIS MOMENT THIS VALUE is MEANINGLESS..
      IF (ISOL.eq.0) SOLACT=0.50000

C--- GEOMETRY CONVERSION ---
C--- Invert Z -> -Z
C--- rotate 40.583 degree in X-Y plane
*-ak- z-> -z will be done again in fillnt
c      DIR(3)=-1.*ORGDIR(3)
      DIR(3)=ORGDIR(3)
      DIR(1)= -(COS(3.141593*2435./60./180.)*ORGDIR(1)
     $       +SIN(3.141593*2435./60./180.)*ORGDIR(2))
      DIR(2)=-(-SIN(3.141593*2435./60./180.)*ORGDIR(1)
     $       +COS(3.141593*2435./60./180.)*ORGDIR(2))
    

*-ak- use honda-flux for energy regions < 100MeV
      IF (ENEU.le.FONE(1))  THEN
         FHONL=fnhonfx(ENEU,DIR,IPAR)
         FNFLUFX=FHONL
         RETURN
      ENDIF
      
*-ak- use honda-flux * norm. factors for energy regions > 100GeV
*     these factors were found by fit_fluka02.dat/honda96_high.dat
*     at 100 GeV with integrated cosz.
*      FFLXHNRM=0.
*      if (ipar.eq.12)  FFLXHNRM=0.80
*      if (ipar.eq.-12) FFLXHNRM=0.75
*      if (ipar.eq.14)  FFLXHNRM=0.88
*      if (ipar.eq.-14) FFLXHNRM=0.82
*-ak below for fit_fluka02.dat/hkkm2.dat
      FFLXHNRM=0.
      if (ipar.eq.12)  FFLXHNRM=1.00
      if (ipar.eq.-12) FFLXHNRM=1.00
      if (ipar.eq.14)  FFLXHNRM=1.00
      if (ipar.eq.-14) FFLXHNRM=1.00
      IF (ENEU.gt.FONE(61)) THEN
         FFLUH=fnhonfx(ENEU,DIR,IPAR)*FFLXHNRM
         FNFLUFX=FFLUH
         RETURN
      ENDIF

      FFLU=FNFDFX(ENEU,DIR,IPAR)
      FNFLUFX=FFLU
      RETURN

      END

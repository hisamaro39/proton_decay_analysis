************************************************************************
*     ---------------------------------
      FUNCTION FNBESFX(ENEU,ORGDIR,IPAR,INDEX)
*     ---------------------------------
C     
*     (Purpose)
C        Honda Differential Flux base on BESS experiment
C
C     (Input)
C       ENEU   : Energy of Neutrino(in GeV)
C       ORGDIR : Direction of the neutrino
C       IPAR   : Flavor of neutrino
C       Index  :       
C                   1 : Flux = Energy^{-2.71}
C                   2 : Flux = Energy^{-2.74}
C
*     (Output)
C       Flux  (m^{-2}sr^{-1}sec^{-1}GeV^{-1})
C
C     (Creation Date and Author)
C      1998.07.    made by J.Kameda based on fnhonfx.F
C                                             
************************************************************************

C
      IMPLICIT NONE
#include"bessflx.h"
#include"nesolact.h"
      
      REAL*4    FNBESFX
      REAL*4    ENEU
      REAL*4    ORGDIR(3)
      INTEGER*4 IPAR
      INTEGER*4 INDEX

      REAL*4    FBESL,FBESM

      REAL*4    FNBSFXL,FNBSFXM,FNBSFXH
      EXTERNAL  FNBSFXL,FNBSFXM,FNBSFXH

      REAL*4    DIR(3)

      INTEGER*4 ISTART
      DATA      ISTART/0/

      IF (ISTART.eq.0) THEN
         ISTART=1
         CALL NERDBES
      ENDIF

C---CHECK COMMON VAR(ISOL == 0 <--> SOLACT = 0.5)
C---AT THIS MOMENT THIS VALUE is MEANINGLESS..
      IF (ISOL.eq.0) SOLACT=0.50000

C---CHECK INDEX
      IF ((INDEX.ne.1).and.(INDEX.ne.2)) then
         write(*,*) "*** ERROR *** Stop executing ***"
         write(*,*) "***       in FNBESFX.F       ***"
         write(*,*) "   Index = ",index
         write(*,*) "   Index should be  1 or 2      "
         STOP
      ENDIF

C--- GEOMETRY CONVERSION ---
C--- Invert Z -> -Z
C--- rotate 40.583 degree in X-Y plane

      DIR(3)=-1.*ORGDIR(3)
      DIR(1)= -(COS(3.141593*2435./60./180.)*ORGDIR(1)
     $       +SIN(3.141593*2435./60./180.)*ORGDIR(2))
      DIR(2)=-(-SIN(3.141593*2435./60./180.)*ORGDIR(1)
     $       +COS(3.141593*2435./60./180.)*ORGDIR(2))
    

      IF (ENEU.lt.BESEM(1))  THEN
         FBESL=FNBSFXL(ENEU,DIR,IPAR)
         FNBESFX=FBESL
         RETURN
      ENDIF

      FBESM=FNBSFXM(ENEU,DIR,IPAR,INDEX)
      FNBESFX=FBESM
      RETURN

      END

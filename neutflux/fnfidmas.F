************************************************************************
*     ---------------------------------
      REAL FUNCTION FNFIDMAS(DWALL)
*     ---------------------------------
*
*     (Purpose)
C        Calculate volume.
*
*     (Input)
C        DWALL : Distance from wall
*
*     (Output)
C        Mass
C
C     (Creation Date and Author)
C        1997.??.?? ; Y.Hayato
C        2007.11.16 ; T.Tanaka
C                     add upmu mode
C        2007.12.05 ; G.Mitsuka
C                     FC/PC and Upmu-water is unified(full water volume)
************************************************************************



      IMPLICIT NONE
#include "necardap.h"
      REAL*4 DWALL
      REAL*4 RAD,HEIGHT,VOL

#include "geotnk.h"
#ifdef UPMU_RCK
#include "upmu.h"
      REAL*4 radrck
      WRITE(*,*) "**** UPMU ROCK ******"
      radrck=DWALL
      VOL=4./3.*PI*(radrck**3)
     &     -PI*(rtktk**2)*hitktk
      if(vol.le.0.) then
        print *,'***Error in fnfidmas; volume=',vol
        stop
      endif
#else
      if (NEWATERVOL.eq.0) then
*     Simulated volume is the entire water tank
         vol=3.1415926535*rtktk*rtktk*hitktk         
         WRITE(*,*) "**** FULL WATER VOLUME ******",vol," cm^3"
      elseif (NEWATERVOL.eq.1) then
         RAD=RINTK - DWALL
         IF (RAD.le.0) THEN
            FNFIDMAS=0.
            RETURN
         ENDIF
         
         HEIGHT=HIINTK-DWALL*2
         IF (HEIGHT.le.0) THEN
            FNFIDMAS=0.
            RETURN
         ENDIF
         VOL=3.1415926535*(RAD**2)*HEIGHT
c              (cm^3) (kg)  (ton) (k)
         WRITE(*,*) "**** ID WATER VOLUME ******",vol," cm^3, DWALL ",DWALL
      elseif (NEWATERVOL.eq.2) then
         vol=3.1415926535*rtktk*rtktk*hitktk
     &        -3.1415926535*rintk*rintk*hiintk
         WRITE(*,*) "**** OD WATER VOLUME ******",vol," cm^3"
      endif
#endif
      FNFIDMAS=VOL   /1000./1000./1000.

      RETURN
      
      END
      

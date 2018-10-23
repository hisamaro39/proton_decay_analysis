      PROGRAM MAIN

      IMPLICIT NONE

#include "hondaflx.h"
#include "nesolact.h"
      INTEGER*4  KIND(4)
      DATA KIND/12,-12,14,-14/
      REAL*4     FLUX(4)
      REAL*4     DIR(3)

      INTEGER*4  IE,I,J,K,ICOS
      REAL*4     ENEU,PHI
      REAL*4     FNHONFX
      EXTERNAL   FNHONFX

      REAL FACTOR

      write(*,*) "FACTOR="
      READ(*,*) FACTOR


      ISOL=1
      SOLACT=0.
      DO 10 K=1,4
         FLUX(K)=0
 10   CONTINUE
         
C      DO 100 IE=06,140
C         ENEU = (10**(REAL((IE/20)-3)))
C         ENEU = ENEU*(10**(REAL(IE-(IE/20)*20)*0.05))
       DO 110 I=1,20
          DO 120 J=1,12
             PHI=(3.141593/6.)*(REAL(J)-0.0)
C             DO 100 IE=41,90
             DO 100 IE=41,90
C             IE=70
                ENEU = (10**(REAL((IE/20)-3)))
                ENEU = ENEU*(10**(REAL(IE-(IE/20)*20)*0.05))*FACTOR
               
C               DIR(3)=-1.05+(REAL(I)*0.1)
               DIR(3)=(REAL(I)*0.1)-0.1000001
               IF (DIR(3).gt.1.0) THEN
                  DIR(3)=0.9-DIR(3)
               ENDIF
               DIR(1)=COS(PHI)*SQRT(1.-DIR(3)**2)
               DIR(2)=SIN(PHI)*SQRT(1.-DIR(3)**2)

               IF (I.gt.11) THEN
                  ICOS=21-I
               ELSE
                  ICOS=I+9
               ENDIF
               
               DO 130 K=1,4
                  FLUX(K)=FNHONFX(ENEU,DIR,KIND(K))
 130           CONTINUE
       IF ((abs((FLUX(3)-HFLXM(IE-40,ICOS,J,3,1))/FLUX(3)).gt.1.E-2)
     $.or.(abs((FLUX(4)-HFLXM(IE-40,ICOS,J,4,1))/FLUX(4).gt.1.E-2))
     $.or.(abs((FLUX(1)-HFLXM(IE-40,ICOS,J,1,1))/FLUX(1).gt.1.E-2))
     $.or.(abs((FLUX(2)-HFLXM(IE-40,ICOS,J,2,1))/FLUX(2).gt.1.E-2)))
     $              THEN
           write(*,'(E11.3,E11.3,E11.3,E11.3,E11.3,E11.3,E11.3)') 
     $          ENEU,DIR(3),PHI,
     $          (FLUX(3)-HFLXM(IE-40,ICOS,J,3,1))/FLUX(3)*100.,
     $          (FLUX(4)-HFLXM(IE-40,ICOS,J,4,1))/FLUX(4)*100.,
     $          (FLUX(1)-HFLXM(IE-40,ICOS,J,1,1))/FLUX(1)*100.,
     $          (FLUX(2)-HFLXM(IE-40,ICOS,J,2,1))/FLUX(2)*100.
           write(999,'(E11.3,E11.3,E11.3,E11.3,E11.3,E11.3,E11.3)') 
     $          ENEU,DIR(3),PHI,FLUX(3),FLUX(4),FLUX(1),FLUX(2)
           write(999,'(E11.3,E11.3,E11.3,E11.3,E11.3,E11.3,E11.3)') 
     $          ENEU,DIR(3),PHI,
     $          HFLXM(IE-40,ICOS,J,3,1),
     $          HFLXM(IE-40,ICOS,J,4,1),
     $          HFLXM(IE-40,ICOS,J,1,1),
     $          HFLXM(IE-40,ICOS,J,2,1)

        ENDIF

 100        CONTINUE
 120     CONTINUE
 110  CONTINUE
         DO 140 K=1,4
            FLUX(K) = FLUX(K) / (4.*3.141593 )
 140     CONTINUE
      END

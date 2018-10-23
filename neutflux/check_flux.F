      PROGRAM MAIN

      IMPLICIT NONE

      INTEGER*4  KIND(4)
      DATA KIND/12,-12,14,-14/
      REAL*4     FLUX(4)
      REAL*4     DIR(3)

      INTEGER*4  IE,I,J,K
      REAL*4     ENEU,PHI
      REAL*4     FNHONFX
      EXTERNAL   FNHONFX

      DO 10 K=1,4
         FLUX(K)=0
 10   CONTINUE
         
      DO 100 IE=06,140
C         write(*,*) "(",IE,"/20)-3=",IE/20-3
C         write(*,*) 
C     $   "REAL(",IE,"-(",IE,"/20)*20=",REAL(IE-(IE/20)*20)*0.05
         ENEU = (10**(REAL((IE/20)-3)))
C         write(*,*) "ENEU=",ENEU
         ENEU = ENEU*(10**(REAL(IE-(IE/20)*20)*0.05))*1.0
C         write(*,*) "ENEU=",ENEU
         DO 110 I=1,20
            DO 120 J=1,12
               PHI=(3.141593/6.)*(REAL(J)-0.5)
               
               DIR(3)=-1.05+(REAL(I)*0.1)
               DIR(1)=COS(PHI)*SQRT(1.-DIR(3)**2)
               DIR(2)=SIN(PHI)*SQRT(1.-DIR(3)**2)

               DO 130 K=1,4
                  FLUX(K)=FLUX(K)+FNHONFX(ENEU,DIR,KIND(K))
 130           CONTINUE
 120        CONTINUE
 110     CONTINUE
         DO 140 K=1,4
            FLUX(K) = FLUX(K) / (4.*3.141593 )
 140     CONTINUE
         write(*,*) ENEU,FLUX(1),FLUX(2),FLUX(3),FLUX(4)
 100  CONTINUE
      END

************************************************************************
*     ---------------------------------
      SUBROUTINE NERNNUDH_1D(ENEU, DIR , IPAR)
*     ---------------------------------
C     
C     (Purpose)
C        Give random  direction of neutrino  by using Honda's flux.
C      
C     (Input)
C        ENEU   : Incident neutrino energy
C        IPAR   : particle code      
C
C     (Output)    
C        DIR    : Direction
C     (Creation Date and Author)
C       1996.??.?? ; Y.Hayato
************************************************************************


      REAL*4    ENEU
      REAL*4    DIR(3)
      INTEGER*4 IPAR

      INTEGER*4 I

      REAL*4    FLXMAP(18,40),FLXMAX,DUM,FLXTMP,DUMRNDM
      
      DO 10 I=1,3
         DIR(I)=0.E0
 10   CONTINUE

      FLXMAX=0.E0
      
      DO 20 I=1,40
         DO 30 J=1,18
            DIR(3)=-1.+REAL(I)*0.05E0
            DIR(1)=COS(REAL(J)*(3.141593/9.))*(sqrt(1.-DIR(3)**2))
            DIR(2)=SIN(REAL(J)*(3.141593/9.))*(sqrt(1.-DIR(3)**2))
            FLXMAP(J,I)=FNHONFX_1D(ENEU,DIR,IPAR)
            IF (FLXMAX.le.FLXMAP(J,I)) FLXMAX=FLXMAP(J,I)
 30      CONTINUE
 20   CONTINUE
      FLXMAX = FLXMAX * 2.

 40   CONTINUE
         DIR(3)=-1.+RLU(DUM)*2.
         PHI  =2.*RLU(DUM)*3.141593
         DIR(1)=COS(PHI)*sqrt(1.-DIR(3)**2)
         DIR(2)=SIN(PHI)*sqrt(1.-DIR(3)**2)
C         DIR(3)=COS(THETA)

         FLXTMP = FNHONFX_1D(ENEU,DIR,IPAR)

         DUMRNDM= RLU(DUM)*FLXMAX
      
      IF (FLXTMP.lt.DUMRNDM) GOTO 40

      RETURN
      
 999  STOP
      END
      

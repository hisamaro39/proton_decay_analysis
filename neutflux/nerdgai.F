************************************************************************
*     ---------------------------------
      SUBROUTINE NERDGAI
*     ---------------------------------
C     
*     (Purpose)
C        PROGRAM to read GAISSER FLUX('96)
C
C     (Input)
C        DATA FILE:gaisser96.dat
C           
C     (Creation Date and Author)
C       1996.??.?? ; Y.Hayato
************************************************************************

      IMPLICIT NONE
#include "gaisser.h"
      
      REAL*4 ELOW(23),EHIGH(40)
      REAL*4 FLEHIGH(2,40,6),FLEBHIGH(2,40,6)
      REAL*4 FLMHIGH(2,40,6),FLMBHIGH(2,40,6)
      
      INTEGER*4 LUNGFLX
      PARAMETER (LUNGFLX=76)

      INTEGER*4 I,J,K,ISOL,IHNDL,KKK
      CHARACTER*80  TMPSTR 
      CHARACTER*256 FILENAME

      INTEGER*4 ILOADED
      SAVE      ILOADED
      DATA      ILOADED/0/

      
C---- CHECK CALLED TWICE?
      IF (ILOADED.ne.0) RETURN
      write(*,*) 'Now LOADING GAISSER-96 FLUX'
      ILOADED = 1
c
c     Rflist style open, kept for backwards compatibility
c
      CALL SKOPENF(LUNGFLX,1,'f',IHNDL)
      IF (IHNDL.LT.0) THEN
c
c        Check environment variable
c
         CALL GETENV("ATMPD_GAISSER",FILENAME)
         IF (FILENAME(1:1) .EQ. ' ') THEN
c
c           Try for the default file with findconsts.
c
            CALL FINDCONSTS(FILENAME,'Gaisser Flux','gaisser96.dat',2)
         END IF
         IF (FILENAME(1:1) .EQ. ' ') THEN
c
c           File still not found. Write error messages below
c           and exit just as before. (IHNDL still < 0)
c
            CONTINUE
         ELSE
            IHNDL = 0
            OPEN(LUNGFLX,FILE=FILENAME,FORM='FORMATTED',STATUS='OLD')
         END IF
      END IF
c
c     Old RFLIST behavior- kept for backwards compatibility...
c
      IF (IHNDL.lt.0) THEN
         IF (IHNDL.eq.-1001) THEN
            write(*,*) "******    ERROR    ******"
            write(*,*) "You did not specify the flux file!"
         ENDIF
         write(*,*) "******    ERROR    ******"
         write(*,*) "Could not open the FLUX file!"
         write(*,*) "LUN=",LUNGFLX,":FILENAME=gaisser96.dat"
         STOP
      ENDIF
         
C--READ LOW ENERGY TABLES
C     K:LOOP OVER SOLAR MIN-MAX (1:MAX,2:MIN)
      DO 30 K=1,2
         IF (K.eq.1) THEN 
            ISOL=2
         ENDIF
         IF (K.eq.2) THEN
            ISOL=1
         ENDIF
C     J:LOOP OVER COSINE
C       1.00(-1.00) , 0.75(-0.75) , 0.50(-0.50) , 0.25(-0.25) ,
C       0.15(-0.15) , 0.05(-0.05)
         DO 20 J=1,6
C     SKIP COMMENT STRING
            read(LUNGFLX,'(A80)') TMPSTR
C            write(*,*) TMPSTR
C     I:LOOP OVER EGERGY
            DO 10 I=1,23
               read(LUNGFLX,*) ELOW(I),
     $              FLGAIS(ISOL,I,J,1),FLGAIS(ISOL,I,J,2),
     $              FLGAIS(ISOL,I,J,3),FLGAIS(ISOL,I,J,4),
     $              FLGAIS(ISOL,I,12+1-J,1),
     $              FLGAIS(ISOL,I,12+1-J,2),
     $              FLGAIS(ISOL,I,12+1-J,3),
     $              FLGAIS(ISOL,I,12+1-J,4)
               EGAIS(I)=10**ELOW(I)
               DO 5 KKK=1,4
                  FLGAIS(ISOL,I,J,KKK)   
     $                 =FLGAIS(ISOL,I,J,KKK)*10000./EGAIS(I)
                  FLGAIS(ISOL,I,13-J,KKK)
     $                 =FLGAIS(ISOL,I,13-J,KKK)*10000./EGAIS(I)
 5             CONTINUE
 10         CONTINUE
 20      CONTINUE
 30   CONTINUE
      
C--READ HIGH ENERGY TABLES(SOL MIN=MAX, COS +=-)
C     J:LOOP OVER COSINE
      DO 50 J=1,6
         read(LUNGFLX,'(A80)') TMPSTR
C         write(*,*) TMPSTR
C     I:LOOP OVER ENERGY
         DO 40 I=1,40
            read(LUNGFLX,*) EHIGH(I),
     $           FLEHIGH(1,I,J),FLEBHIGH(1,I,J),
     $           FLMHIGH(1,I,J),FLMBHIGH(1,I,J)
            EGAIS   (I+23)=EHIGH(I)
            FLGAIS  (1,I+23,J     ,1)=FLEHIGH(1,I,J)*10000./EHIGH(I)
            FLGAIS  (1,I+23,12+1-J,1)=FLEHIGH(1,I,J)*10000./EHIGH(I)
            FLGAIS  (2,I+23,J     ,1)=FLEHIGH(1,I,J)*10000./EHIGH(I)
            FLGAIS  (2,I+23,12+1-J,1)=FLEHIGH(1,I,J)*10000./EHIGH(I)
            
            FLGAIS  (1,I+23,J     ,2)=FLEBHIGH(1,I,J)*10000./EHIGH(I)
            FLGAIS  (1,I+23,12+1-J,2)=FLEBHIGH(1,I,J)*10000./EHIGH(I)
            FLGAIS  (2,I+23,J     ,2)=FLEBHIGH(1,I,J)*10000./EHIGH(I)
            FLGAIS  (2,I+23,12+1-J,2)=FLEBHIGH(1,I,J)*10000./EHIGH(I)

            FLGAIS  (1,I+23,J     ,3)=FLMHIGH(1,I,J)*10000./EHIGH(I)
            FLGAIS  (1,I+23,12+1-J,3)=FLMHIGH(1,I,J)*10000./EHIGH(I)
            FLGAIS  (2,I+23,J     ,3)=FLMHIGH(1,I,J)*10000./EHIGH(I)
            FLGAIS  (2,I+23,12+1-J,3)=FLMHIGH(1,I,J)*10000./EHIGH(I)
            
            FLGAIS  (1,I+23,J     ,4)=FLMBHIGH(1,I,J)*10000./EHIGH(I)
            FLGAIS  (1,I+23,12+1-J,4)=FLMBHIGH(1,I,J)*10000./EHIGH(I)
            FLGAIS  (2,I+23,J     ,4)=FLMBHIGH(1,I,J)*10000./EHIGH(I)
            FLGAIS  (2,I+23,12+1-J,4)=FLMBHIGH(1,I,J)*10000./EHIGH(I)
            
 40      CONTINUE
 50   CONTINUE
      CALL SKCLOSEF(LUNGFLX)

      DO 100 ISOL=1,2
         DO 110 K=1,4
            DO 120 I=1,63
               FLGAIST(ISOL,I,K)=0.
               DO 130 J=1,11
                  FLGAIST(ISOL,I,K)=FLGAIST(ISOL,I,K)+
     $            ( (FLGAIS(ISOL,I,J,K)+FLGAIS(ISOL,I,J,K))
     $             *(COSGAI(J)-COSGAI(J+1))/2.)*2*3.141592
 130           CONTINUE
               FLGAIST(ISOL,I,K)=FLGAIST(ISOL,I,K)/4./3.141593
 120        CONTINUE
 110     CONTINUE
 100  CONTINUE

      RETURN
      END

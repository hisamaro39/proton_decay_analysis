************************************************************************
*     ---------------------------------
      SUBROUTINE NERDBES
*     ---------------------------------
C     
*     (Purpose)
C        READ FLUX CALCULATED BY M.Honda
C           
C     (Input) 
C        FLUX FOR LOW  ENERGY(  2MeV to 112MeV)/temporary KAMLOW
C        FLUX FOR MID  ENERGY(112MeV to  30GeV)
C        FLUX FOR HIGH ENERGY( 30GeV to   1TeV)
C     (Creation Date and Author)
C       1996.??.?? ; Y.Hayato
C
C
************************************************************************

      IMPLICIT NONE
      
#include"bessflx.h"
C      REAL*4 BFLXL(72,20,4,2),BFLXM(50,20,12,4,2),BFLXH(31,20,4),
C     $       BFLX(100,4,2),BFLXA(153,4,2)
C      REAL*4 BESE(100),BESEL(72),BESEM(50),BESEH(31)
C      REAL*4 BESEA(72+50+31)
C      REAL*4 BFLXLNRM(4,2)
C      COMMON /BESFX/BESE(35+50+15) ,BFLX(100,4,2),
C     $              BESEA(72+50+31),BFLXA(153,4,2),
C     $              BESEL(72),BFLXL(72,20,4,2),
C     $              BESEM(50),BFLXM(50,20,12,4,2),
C     $              BESEH(31),BFLXH(31,20,4),
C     $              BFLXLNRM(4,2)


      INTEGER*4 LUNBFLXL
      PARAMETER (LUNBFLXL=87)
      INTEGER*4 LUNBFLX1
      PARAMETER (LUNBFLX1=88)
      INTEGER*4 LUNBFLX2
      PARAMETER (LUNBFLX2=89)

      INTEGER*4 IHNDL

      INTEGER*4 IK,IE,IZ,IA,ICOS,IS,INDEX
      CHARACTER*80 DUMSTR
      REAL*4    DUM(10),BFLXLTMP(4),DUME

      INTEGER*4 ILOADED
      SAVE      ILOADED
      DATA      ILOADED/0/

      
C---- CHECK CALLED TWICE?
      IF (ILOADED.ne.0) RETURN
      write(*,*) "Now LOADING HONDA FLUX BASED ON BESS RESULTS"
      ILOADED = 1

C---  Normalize Factor for L.E.Table
C     Temporary(To adjust old-table and new table)
C     Factor was decided to match total flax

      BFLXLNRM(1,1)=1.00
      BFLXLNRM(2,1)=1.05
      BFLXLNRM(3,1)=1.00
      BFLXLNRM(4,1)=1.00

      BFLXLNRM(1,2)=0.89
      BFLXLNRM(2,2)=0.95
      BFLXLNRM(3,2)=0.89
      BFLXLNRM(4,2)=0.89

      CALL SKOPENF(LUNBFLXL,1,'f',IHNDL)
      IF (IHNDL.lt.0) THEN
         IF (IHNDL.eq.-1001) THEN
            write(*,*) "******    ERROR    ******"
            write(*,*) "You did not specify the flux file!"
         ENDIF
         write(*,*) "******    ERROR    ******"
         write(*,*) "Could not open the FLUX file!"
         write(*,*) "LUN=",LUNBFLXL,":FILENAME=honda96low.dat"
         STOP
      ENDIF

      DO 10 IZ=1,20
         READ(LUNBFLXL,'(A80)') DUMSTR
C         write(*,*) "Now reading ",DUMSTR
         DO 20 IE=1,72
            IF (IZ.le.11) THEN
               ICOS=IZ+9
            ELSE
               ICOS=21-IZ
            ENDIF
            READ(LUNBFLXL,*) BESEL(IE),DUM(1),DUM(2),DUM(3),DUM(4),
     $                              BFLXLTMP(3),
     $                              BFLXLTMP(4),
     $                              BFLXLTMP(1),
     $                              BFLXLTMP(2),
     $                              DUM(5),DUM(6),DUM(7)
         BFLXL(IE,ICOS,1,1)=BFLXLTMP(1)*BFLXLNRM(1,1)
         BFLXL(IE,ICOS,2,1)=BFLXLTMP(2)*BFLXLNRM(2,1)
         BFLXL(IE,ICOS,3,1)=BFLXLTMP(3)*BFLXLNRM(3,1)
         BFLXL(IE,ICOS,4,1)=BFLXLTMP(4)*BFLXLNRM(4,1)

         BFLXL(IE,ICOS,1,2)=BFLXLTMP(1)*BFLXLNRM(1,2)
         BFLXL(IE,ICOS,2,2)=BFLXLTMP(2)*BFLXLNRM(2,2)
         BFLXL(IE,ICOS,3,2)=BFLXLTMP(3)*BFLXLNRM(3,2)
         BFLXL(IE,ICOS,4,2)=BFLXLTMP(4)*BFLXLNRM(4,2)
 20      CONTINUE
C      write(*,*) "IZ=",IZ," is OVER"
 10   CONTINUE

      CALL SKCLOSEF(LUNBFLXL)

      CALL SKOPENF(LUNBFLX1,1,'f',IHNDL)
      IF (IHNDL.lt.0) THEN
         IF (IHNDL.eq.-1001) THEN
            write(*,*) "******    ERROR    ******"
            write(*,*) "You did not specify the flux file!"
         ENDIF
         write(*,*) "******    ERROR    ******"
         write(*,*) "Could not open the FLUX file!"
         write(*,*) "LUN=",LUNBFLX1,":FILENAME=bess71.dat"
         STOP
      ENDIF


         DO 50 IZ=1,20

               READ(LUNBFLX1,'(A80)') DUMSTR

               DO 70 IE=1,81
                  
                  IF (IZ.le.11) THEN
                     ICOS=IZ+9
                  ELSE
                     ICOS=21-IZ
                  ENDIF
                  READ(LUNBFLX1,*) BESEM(IE),
     $                                    BFLXM(IE,ICOS,1,3,1),
     $                                    BFLXM(IE,ICOS,1,4,1),
     $                                    BFLXM(IE,ICOS,1,1,1),
     $                                    BFLXM(IE,ICOS,1,2,1)

            DO 60 IA=2,12
                  BFLXM(IE,ICOS,IA,3,1) = BFLXM(IE,ICOS,1,3,1)
                  BFLXM(IE,ICOS,IA,4,1) = BFLXM(IE,ICOS,1,4,1)
                  BFLXM(IE,ICOS,IA,1,1) = BFLXM(IE,ICOS,1,1,1)
                  BFLXM(IE,ICOS,IA,2,1) = BFLXM(IE,ICOS,1,2,1)
 60         CONTINUE


C                  write(*,'(E11.4,E11.4,E11.4,E11.4,E11.4)')
C     $                  BESEM(IE),
C     $                  BFLXM(IE,ICOS,IA,3,IS),
C     $                  BFLXM(IE,ICOS,IA,4,IS),
C     $                  BFLXM(IE,ICOS,IA,1,IS),
C     $                  BFLXM(IE,ICOS,IA,2,IS)


 70            CONTINUE

 50      CONTINUE

      CALL SKCLOSEF(LUNBFLX1)

      CALL SKOPENF(LUNBFLX2,1,'f',IHNDL)
      IF (IHNDL.lt.0) THEN
         IF (IHNDL.eq.-1001) THEN
            write(*,*) "******    ERROR    ******"
            write(*,*) "You did not specify the flux file!"
         ENDIF
         write(*,*) "******    ERROR    ******"
         write(*,*) "Could not open the FLUX file!"
         write(*,*) "LUN=",LUNBFLX2,":FILENAME=bess71.dat"
         STOP
      ENDIF


         DO 150 IZ=1,20

               READ(LUNBFLX2,'(A80)') DUMSTR

               DO 170 IE=1,81
                  
                  IF (IZ.le.11) THEN
                     ICOS=IZ+9
                  ELSE
                     ICOS=21-IZ
                  ENDIF
                  READ(LUNBFLX2,*) DUME,
     $                                    BFLXM(IE,ICOS,1,3,2),
     $                                    BFLXM(IE,ICOS,1,4,2),
     $                                    BFLXM(IE,ICOS,1,1,2),
     $                                    BFLXM(IE,ICOS,1,2,2)

            DO 160 IA=2,12
                  BFLXM(IE,ICOS,IA,3,2) = BFLXM(IE,ICOS,1,3,2)
                  BFLXM(IE,ICOS,IA,4,2) = BFLXM(IE,ICOS,1,4,2)
                  BFLXM(IE,ICOS,IA,1,2) = BFLXM(IE,ICOS,1,1,2)
                  BFLXM(IE,ICOS,IA,2,2) = BFLXM(IE,ICOS,1,2,2)
 160         CONTINUE


C                  write(*,'(E11.4,E11.4,E11.4,E11.4,E11.4)')
C     $                  BESEM(IE),
C     $                  BFLXM(IE,ICOS,IA,3,IS),
C     $                  BFLXM(IE,ICOS,IA,4,IS),
C     $                  BFLXM(IE,ICOS,IA,1,IS),
C     $                  BFLXM(IE,ICOS,IA,2,IS)


 170            CONTINUE

 150      CONTINUE

      CALL SKCLOSEF(LUNBFLX2)


CC ----------- Integration -------

      DO 190 INDEX = 1,2
      DO 200 IS=1,2
         DO 210 IK=1,4
            DO 220 IE=1,115
               BFLX(IE,IK,IS,INDEX)=0.
               IF (IE.le.34) THEN
                  BESE(IE)=BESEL(IE)
                  DO 230 IZ=1,20
                     BFLX(IE,IK,IS,INDEX)= BFLX(IE,IK,IS,INDEX)
     $                              +BFLXL(IE,IZ,IK,IS)
 230              CONTINUE
           BFLX(IE,IK,IS,INDEX)=BFLX(IE,IK,IS,INDEX)*0.1E0*2*3.141592
               ENDIF

               IF (IE.ge.35) THEN
                  BESE(IE)=BESEM(IE-34)
                  DO 240 IZ=1,20
                     DO 250 IA=1,12
                        BFLX(IE,IK,IS,INDEX)= BFLX(IE,IK,IS,INDEX)
     $                                 +BFLXM(IE-34,IZ,IA,IK,INDEX)
 250                 CONTINUE
 240              CONTINUE
         BFLX(IE,IK,IS,INDEX)=BFLX(IE,IK,IS,INDEX)*0.1E0*(3.141592/6.)
               ENDIF

               BFLX(IE,IK,IS,INDEX)=BFLX(IE,IK,IS,INDEX)/(4.*3.141592)
 220        CONTINUE
 210     CONTINUE
 200  CONTINUE
 190  CONTINUE

      DO 290 INDEX = 1,2
      DO 300 IS=1,2
         DO 310 IK=1,4
            DO 320 IE=1,153
               BFLXA(IE,IK,IS,INDEX)=0.
               IF (IE.le.72) THEN
                  BESEA(IE)=BESEL(IE)
                  DO 330 IZ=1,20
                     BFLXA(IE,IK,IS,INDEX)= BFLXA(IE,IK,IS,INDEX)
     $                               +BFLXL(IE,IZ,IK,IS)
 330              CONTINUE
            BFLXA(IE,IK,IS,INDEX)=BFLXA(IE,IK,IS,INDEX)*0.1E0*2*3.141592
               ENDIF
               IF (IE.ge.73) THEN
                  BESEA(IE)=BESEM(IE-72)
                  DO 340 IZ=1,20
                     DO 350 IA=1,12
                        BFLXA(IE,IK,IS,INDEX)= BFLXA(IE,IK,IS,INDEX)
     $                                 +BFLXM(IE-72,IZ,IA,IK,INDEX)
 350                 CONTINUE
 340              CONTINUE
         BFLXA(IE,IK,IS,INDEX)=BFLXA(IE,IK,IS,INDEX)*0.1E0*(3.141592/6.)
               ENDIF

 320        CONTINUE
 310     CONTINUE
 300  CONTINUE
 290  CONTINUE

      END


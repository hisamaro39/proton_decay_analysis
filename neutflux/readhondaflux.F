      SUBROUTINE RDHONFX

      IMPLICIT NONE
      
C     READ FLUX CALCULATED BY M.Honda
C           
C     FLUX FOR LOW  ENERGY(  2MeV to 112MeV)/temporary KAMLOW
C     FLUX FOR MID  ENERGY(112MeV to  30GeV)/
C     FLUX FOR HIGH ENERGY( 30GeV to   1TeV)
C      
C
      
#include"hondaflx.h"
C      REAL*4 HFLXL(72,20,4),HFLXM(50,20,12,4,2),HFLXH(31,20,4),
C     $       HFLX(100,4,2),HFLXA(153,4,2)
C      REAL*4 HONE(100),HONEL(72),HONEM(50),HONEH(31)
C      REAL*4 HONEA(72+50+31)
C      COMMON /HONFX/HONE(35+50+15) ,HFLX(100,4,2),
C     $              HONEA(72+50+31),HFLXA(153,4,2),
C     $              HONEL(72)      ,HFLXL(72,20,4),
C     $              HONEM(50)      ,HFLXM(50,20,12,4,2),
C     $              HONEH(31)      ,HFLXH(31,20,4)


      INTEGER*4 LUNHFLXL
      PARAMETER (LUNHFLXL=77)
      INTEGER*4 LUNHFLXM
      PARAMETER (LUNHFLXM=78)
      INTEGER*4 LUNHFLXH
      PARAMETER (LUNHFLXH=79)

      INTEGER*4 IHNDL

      INTEGER*4 IK,IE,IZ,IA,ICOS,IS
      CHARACTER*80 DUMSTR
      REAL*4    DUM(10)

      REAL*4 DUMTMP,HONEH(31),HFLXH(31,20,4)
C     CALL SKOPENF(LUNHFLXL,0,'C',IHNDL)
C     CALL SKOPENF(LUNHFLXM,0,'C',IHNDL)
C     CALL SKOPENF(LUNHFLXH,0,'C',IHNDL)

      DO 10 IZ=1,20
         READ(LUNHFLXL,'(A80)') DUMSTR
C         write(*,*) "Now reading ",DUMSTR
         DO 20 IE=1,72
            IF (IZ.le.11) THEN
               ICOS=IZ+9
            ELSE
               ICOS=21-IZ
            ENDIF
            READ(LUNHFLXL,*) HONEL(IE),DUM(1),DUM(2),DUM(3),DUM(4),
     $                              HFLXL(IE,ICOS,3,dumtmp),
     $                              HFLXL(IE,ICOS,4,dumtmp),
     $                              HFLXL(IE,ICOS,1,dumtmp),
     $                              HFLXL(IE,ICOS,2,dumtmp),
     $                              DUM(5),DUM(6),DUM(7)
C         HFLXL(IE,ICOS,1)=HFLXL(IE,ICOS,1)/HONEL(IE)
C         HFLXL(IE,ICOS,2)=HFLXL(IE,ICOS,2)/HONEL(IE)
C         HFLXL(IE,ICOS,3)=HFLXL(IE,ICOS,3)/HONEL(IE)
C         HFLXL(IE,ICOS,4)=HFLXL(IE,ICOS,4)/HONEL(IE)
 20      CONTINUE
C      write(*,*) "IZ=",IZ," is OVER"
 10   CONTINUE

      DO 40 IS=1,2
         DO 50 IZ=1,20
            DO 60 IA=1,12
               READ(LUNHFLXM,'(A80)') DUMSTR
               write(*,*) "Now reading ",DUMSTR
               DO 70 IE=1,50
                  IF (IZ.le.11) THEN
                     ICOS=IZ+9
                  ELSE
                     ICOS=21-IZ
                  ENDIF
                  READ(LUNHFLXM,*) HONEM(IE),HFLXM(IE,ICOS,IA,3,IS),
     $                                    HFLXM(IE,ICOS,IA,4,IS),
     $                                    HFLXM(IE,ICOS,IA,1,IS),
     $                                    HFLXM(IE,ICOS,IA,2,IS)
                  
 70            CONTINUE
 60         CONTINUE
 50      CONTINUE
 40   CONTINUE

      DO 80 IK=1,4
         READ(LUNHFLXH,'(A80)') DUMSTR
         DO 90 IE=1,31
            READ(LUNHFLXH,*)          HONEH(IE),HFLXH(IE,11,IK),
     $           HFLXH(IE,12,IK),HFLXH(IE,13,IK),HFLXH(IE,14,IK),
     $           HFLXH(IE,15,IK),HFLXH(IE,16,IK),HFLXH(IE,17,IK),
     $           HFLXH(IE,18,IK),HFLXH(IE,19,IK),HFLXH(IE,20,IK)
            DO 100 IZ=1,10
               HFLXH(IE,IZ+10,IK)=HFLXH(IE,IZ+10,IK)/(HONEH(IE)**3)
               HFLXH(IE,11-IZ,IK)=HFLXH(IE,IZ+10,IK)
 100        CONTINUE
 90      CONTINUE
 80   CONTINUE

      DO 200 IS=1,2
         DO 210 IK=1,4
            DO 220 IE=1,100
               HFLX(IE,IK,IS)=0.
               IF (IE.le.35) THEN
                  HONE(IE)=HONEL(IE)
                  DO 230 IZ=1,20
                     HFLX(IE,IK,IS)=HFLX(IE,IK,IS)+HFLXL(IE,IZ,IK,dumtmp)
 230              CONTINUE
                  HFLX(IE,IK,IS)=HFLX(IE,IK,IS)*0.1E0*2*3.141592
               ENDIF
               IF (IE.ge.36.and.IE.le.85) THEN
                  HONE(IE)=HONEM(IE-35)
                  DO 240 IZ=1,20
                     DO 250 IA=1,12
                        HFLX(IE,IK,IS)= HFLX(IE,IK,IS)
     $                                 +HFLXM(IE-35,IZ,IA,IK,IS)
 250                 CONTINUE
 240              CONTINUE
                  HFLX(IE,IK,IS)=HFLX(IE,IK,IS)*0.1E0*(3.141592/6.)
               ENDIF
               IF (IE.ge.86) THEN
                  HONE(IE)=HONEH(IE-69)
                  DO 260 IZ=1,20
                     HFLX(IE,IK,IS)=  HFLX(IE,IK,IS)
     $                              +HFLXH(IE-69,IZ,IK)
 260              CONTINUE
                  HFLX(IE,IK,IS)=HFLX(IE,IK,IS)*0.1E0*2*3.141592
               ENDIF
 220        CONTINUE
 210     CONTINUE
 200  CONTINUE

      DO 300 IS=1,2
         DO 310 IK=1,4
            DO 320 IE=1,153
               HFLXA(IE,IK,IS)=0.
               IF (IE.le.72) THEN
                  HONEA(IE)=HONEL(IE)
                  DO 330 IZ=1,20
                     HFLXA(IE,IK,IS)=HFLXA(IE,IK,IS)+HFLXL(IE,IZ,IK,dumtmp)
 330              CONTINUE
                  HFLXA(IE,IK,IS)=HFLXA(IE,IK,IS)*0.1E0*2*3.141592
               ENDIF
               IF (IE.ge.73.and.IE.le.122) THEN
                  HONEA(IE)=HONEM(IE-72)
                  DO 340 IZ=1,20
                     DO 350 IA=1,12
                        HFLXA(IE,IK,IS)= HFLXA(IE,IK,IS)
     $                                 +HFLXM(IE-72,IZ,IA,IK,IS)
 350                 CONTINUE
 340              CONTINUE
                  HFLXA(IE,IK,IS)=HFLXA(IE,IK,IS)*0.1E0*(3.141592/6.)
               ENDIF
               IF (IE.ge.123) THEN
                  HONEA(IE)=HONEH(IE-122)
                  DO 360 IZ=1,20
                     HFLXA(IE,IK,IS)=  HFLXA(IE,IK,IS)
     $                              +HFLXH(IE-122,IZ,IK)
 360              CONTINUE
                  HFLXA(IE,IK,IS)=HFLXA(IE,IK,IS)*0.1E0*2*3.141592
               ENDIF
 320        CONTINUE
 310     CONTINUE
 300  CONTINUE

      DO 900 IS=1,2
         DO 910 IE=1,50
            DO 920 ICOS=1,20
               DO 930 IA=1,12
                  write(*,'(F11.4,F11.4,F11.4,F11.4,F11.4)')
     $                 HONEM(IE),
     $                 HFLXM(IE,ICOS,IA,3,IS),
     $                 HFLXM(IE,ICOS,IA,4,IS),
     $                 HFLXM(IE,ICOS,IA,1,IS),
     $                 HFLXM(IE,ICOS,IA,2,IS)
 930           CONTINUE
 920        CONTINUE
 910     CONTINUE
 900  CONTINUE
      STOP
      END

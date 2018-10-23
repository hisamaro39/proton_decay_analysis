      PROGRAM READHONDAFLUX
C      SUBROUTINE READHONDAFLUX

      IMPLICIT NONE
      
C     READ FLUX CALCULATED BY M.Honda
C           
C     FLUX FOR LOW  ENERGY(  2MeV to 112MeV)/temporary KAMLOW
C     FLUX FOR MID  ENERGY(112MeV to  30GeV)/
C     FLUX FOR HIGH ENERGY( 30GeV to   1TeV)
C      
C
      
      REAL*4 HFLXL(72,20,4),HFLXM(50,20,12,4,2),HFLXH(31,20,4),
     $       HFLX(100,4,2),HFLXA(153,4,2)
      REAL*4 HONE(100),HONEL(72),HONEM(50),HONEH(31)
      REAL*4 HONEA(72+50+31)
      COMMON /HONFX/HONE,HFLX,HFLXA,
     $              HONEL,HFLXL,
     $              HONEM,HFLXM,
     $              HONEH,HFLXH


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
      
C     CALL SKOPENF(LUNHFLXL,0,'C',IHNDL)
C     CALL SKOPENF(LUNHFLXM,0,'C',IHNDL)
C     CALL SKOPENF(LUNHFLXH,0,'C',IHNDL)

      DO 10 IZ=1,20
         READ(LUNHFLXL,'(A80)') DUMSTR
         write(*,*) "Now reading ",DUMSTR
         DO 20 IE=1,72
            IF (IZ.le.11) THEN
               ICOS=IZ+9
            ELSE
               ICOS=21-IZ
            ENDIF
            READ(LUNHFLXL,*) HONEL(IE),DUM(1),DUM(2),DUM(3),DUM(4),
     $                              HFLXL(IE,ICOS,3),
     $                              HFLXL(IE,ICOS,4),
     $                              HFLXL(IE,ICOS,1),
     $                              HFLXL(IE,ICOS,2),
     $                              DUM(5),DUM(6),DUM(7)
C         HFLXL(IE,ICOS,1)=HFLXL(IE,ICOS,1)/HONEL(IE)
C         HFLXL(IE,ICOS,2)=HFLXL(IE,ICOS,2)/HONEL(IE)
C         HFLXL(IE,ICOS,3)=HFLXL(IE,ICOS,3)/HONEL(IE)
C         HFLXL(IE,ICOS,4)=HFLXL(IE,ICOS,4)/HONEL(IE)
 20      CONTINUE
      write(*,*) "IZ=",IZ," is OVER"
 10   CONTINUE

      write(*,*) "LOW ENERGY is OVER"

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
                  write(100,
     $                 '(E11.3,E11.3,E11.3,E11.3,E11.3,E11.3,E11.3)')
     $                          HONEM(IE),REAL(IZ-11)*0.1+0.5,
     $                                    (REAL(IA)-0.5)*(3.141593/6.),
     $                                    HFLXM(IE,ICOS,IA,3,IS),
     $                                    HFLXM(IE,ICOS,IA,4,IS),
     $                                    HFLXM(IE,ICOS,IA,1,IS),
     $                                    HFLXM(IE,ICOS,IA,2,IS)

 70            CONTINUE
 60         CONTINUE
 50      CONTINUE
 40   CONTINUE

      DO 80 IK=1,4
         READ(LUNHFLXH,'(A80)') DUMSTR
         write(*,*) "Now reading ",DUMSTR
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
                     HFLX(IE,IK,IS)=HFLX(IE,IK,IS)+HFLXL(IE,IZ,IK)
 230              CONTINUE
                  HFLX(IE,IK,IS)=HFLX(IE,IK,IS)*0.1E0*2*3.141593
               ENDIF
               IF (IE.ge.36.and.IE.le.85) THEN
                  HONE(IE)=HONEM(IE-35)
                  DO 240 IZ=1,20
                     DO 250 IA=1,12
                        HFLX(IE,IK,IS)= HFLX(IE,IK,IS)
     $                                 +HFLXM(IE-35,IZ,IA,IK,IS)
                     IF (HFLXM(IE-35,IZ,IA,IK,IS).le.0.00001) THEN
                        write(*,*) "HFLXM(",IE-35,",",IZ,",",IA,",",
     $                             IK,",",IS,")=",
     $                  HFLXM(IE-35,IZ,IA,IK,IS)
                     ENDIF
 250                 CONTINUE
 240              CONTINUE
                  HFLX(IE,IK,IS)=HFLX(IE,IK,IS)*0.1E0*(3.141593/6.)
               ENDIF
               IF (IE.ge.86) THEN
                  HONE(IE)=HONEH(IE-69)
                  DO 260 IZ=1,20
                     HFLX(IE,IK,IS)=  HFLX(IE,IK,IS)
     $                              +HFLXH(IE-69,IZ,IK)
 260              CONTINUE
                  HFLX(IE,IK,IS)=HFLX(IE,IK,IS)*0.1E0*2*3.141593
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
                     HFLXA(IE,IK,IS)=HFLXA(IE,IK,IS)+HFLXL(IE,IZ,IK)
 330              CONTINUE
                  HFLXA(IE,IK,IS)=HFLXA(IE,IK,IS)*0.1E0*2*3.141593
               ENDIF
               IF (IE.ge.73.and.IE.le.122) THEN
                  HONEA(IE)=HONEM(IE-72)
                  DO 340 IZ=1,20
                     DO 350 IA=1,12
                        HFLXA(IE,IK,IS)= HFLXA(IE,IK,IS)
     $                                 +HFLXM(IE-72,IZ,IA,IK,IS)
                     IF (HFLXM(IE-72,IZ,IA,IK,IS).le.0.00001) THEN
                        write(*,*) "HFLXM(",IE-72,",",IZ,",",IA,",",
     $                             IK,",",IS,")=",
     $                  HFLXM(IE-72,IZ,IA,IK,IS)
                     ENDIF
 350                 CONTINUE
 340              CONTINUE
                  HFLXA(IE,IK,IS)=HFLXA(IE,IK,IS)*0.1E0*(3.141593/6.)
               ENDIF
               IF (IE.ge.123) THEN
                  HONEA(IE)=HONEH(IE-122)
                  DO 360 IZ=1,20
                     HFLXA(IE,IK,IS)=  HFLXA(IE,IK,IS)
     $                              +HFLXH(IE-122,IZ,IK)
 360              CONTINUE
                  HFLXA(IE,IK,IS)=HFLXA(IE,IK,IS)*0.1E0*2*3.141593
               ENDIF
 320        CONTINUE
 310     CONTINUE
 300  CONTINUE

      GOTO  501
      
      DO 400 IE=1,100
         write(500,*) HONE(IE),HFLX(IE,1,1),HFLX(IE,1,2),
     $                         HFLX(IE,2,1),HFLX(IE,2,2),
     $                         HFLX(IE,3,1),HFLX(IE,3,2),
     $                         HFLX(IE,4,1),HFLX(IE,4,2)
         write(502,*) HONE(IE),HFLX(IE,1,1)*HONE(IE)**3/4./3.14,
     $                         HFLX(IE,1,2)*HONE(IE)**3/4./3.14,
     $                         HFLX(IE,2,1)*HONE(IE)**3/4./3.14,
     $                         HFLX(IE,2,2)*HONE(IE)**3/4./3.14,
     $                         HFLX(IE,3,1)*HONE(IE)**3/4./3.14,
     $                         HFLX(IE,3,2)*HONE(IE)**3/4./3.14,
     $                         HFLX(IE,4,1)*HONE(IE)**3/4./3.14,
     $                         HFLX(IE,4,2)*HONE(IE)**3/4./3.13
 400  CONTINUE

      DO 500 IE=1,153
         write(600,*) HONEA(IE),HFLXA(IE,1,1),HFLXA(IE,1,2),
     $                         HFLXA(IE,2,1),HFLXA(IE,2,2),
     $                         HFLXA(IE,3,1),HFLXA(IE,3,2),
     $                         HFLXA(IE,4,1),HFLXA(IE,4,2)
         write(601,*) HONEA(IE),
     $                HFLXA(IE,1,1)*HONEA(IE)**3/4./3.14,
     $                HFLXA(IE,1,2)*HONEA(IE)**3/4./3.14,
     $                HFLXA(IE,2,1)*HONEA(IE)**3/4./3.14,
     $                HFLXA(IE,2,2)*HONEA(IE)**3/4./3.14,
     $                HFLXA(IE,3,1)*HONEA(IE)**3/4./3.14,
     $                HFLXA(IE,3,2)*HONEA(IE)**3/4./3.14,
     $                HFLXA(IE,4,1)*HONEA(IE)**3/4./3.14,
     $                HFLXA(IE,4,2)*HONEA(IE)**3/4./3.14
         write(602,*) HONEA(IE),HFLXA(IE,1,1)/4./3.14,
     $                         HFLXA(IE,1,2)/4./3.14,
     $                         HFLXA(IE,2,1)/4./3.14,
     $                         HFLXA(IE,2,2)/4./3.14,
     $                         HFLXA(IE,3,1)/4./3.14,
     $                         HFLXA(IE,3,2)/4./3.14,
     $                         HFLXA(IE,4,1)/4./3.14,
     $                         HFLXA(IE,4,2)/4./3.13
 500  CONTINUE
 501  CONTINUE

      END

************************************************************************
*     ---------------------------------
      SUBROUTINE NERDHON_1D
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
************************************************************************

      IMPLICIT NONE
      
#include"hondaflx_1d.h"
C      REAL*4 HFLXL(72,20,4,2),HFLXM(50,20,12,4,2),HFLXH(31,20,4),
C     $       HFLX(100,4,2),HFLXA(153,4,2)
C      REAL*4 HONE(100),HONEL(72),HONEM(50),HONEH(31)
C      REAL*4 HONEA(72+50+31)
C      REAL*4 HFLXLNRM(4,2)
C      COMMON /HONFX/HONE(35+50+15) ,HFLX(100,4,2),
C     $              HONEA(72+50+31),HFLXA(153,4,2),
C     $              HONEL(72),HFLXL(72,20,4,2),
C     $              HONEM(50),HFLXM(50,20,12,4,2),
C     $              HONEH(31),HFLXH(31,20,4),
C     $              HFLXLNRM(4,2)


      INTEGER*4 LUNHFLXL
      PARAMETER (LUNHFLXL=71)
      INTEGER*4 LUNHFLXM
      PARAMETER (LUNHFLXM=72)
      INTEGER*4 LUNHFLXH
      PARAMETER (LUNHFLXH=73)

      INTEGER*4 IHNDL

      INTEGER*4 IK,IE,IZ,IA,ICOS,IS
      CHARACTER*80  DUMSTR 
      CHARACTER*256 FILENAME
      REAL*4    DUM(10),HFLXLTMP(4)

      INTEGER*4 ILOADED
      SAVE      ILOADED
      DATA      ILOADED/0/

      INTEGER*4 IUSE_OPEN
      
C---- CHECK CALLED TWICE?
      IF (ILOADED.ne.0) RETURN
      write(*,*) "Now LOADING HONDA(1D) FLUX"
      ILOADED = 1

C---  Normalize Factor for L.E.Table
C     Temporary(To adjust old-table and new table)
C     Factor was decided to match total flax

      HFLXLNRM_1D(1,1)=1.00
      HFLXLNRM_1D(2,1)=1.05
      HFLXLNRM_1D(3,1)=1.00
      HFLXLNRM_1D(4,1)=1.00

      HFLXLNRM_1D(1,2)=0.89
      HFLXLNRM_1D(2,2)=0.95
      HFLXLNRM_1D(3,2)=0.89
      HFLXLNRM_1D(4,2)=0.89
c
c     Rflist style open, kept for backwards compatibility
c
      IUSE_OPEN=0
      CALL SKOPENF(LUNHFLXL,1,'f',IHNDL)
      IF (IHNDL.LT.0) THEN
c
c        Check environment variable
c
         CALL GETENV("ATMPD_HONDALOW",FILENAME)
         IF (FILENAME(1:1) .EQ. ' ') THEN
c
c           Try for the default file with findconsts.
c
            CALL FINDCONSTS(FILENAME,'Honda Flux Low','honda96low.dat',2)
         END IF
         IF (FILENAME(1:1) .EQ. ' ') THEN
c
c           File still not found. Write error messages below
c           and exit just as before. (IHNDL still < 0)
c
            CONTINUE
         ELSE
            IHNDL = 0
            OPEN(LUNHFLXL,FILE=FILENAME,FORM='FORMATTED',STATUS='OLD')
            IUSE_OPEN=1
         END IF
      END IF
c
c     Old RFLIST behavior- kept for backwards compatibility...
c
      IF (IHNDL.LT.0) THEN
         IF (IHNDL.eq.-1001) THEN
            write(*,*) "******    ERROR    ******"
            write(*,*) "You did not specify the flux file!"
         ENDIF
         write(*,*) "******    ERROR    ******"
         write(*,*) "Could not open the FLUX file!"
         write(*,*) "LUN=",LUNHFLXL,":FILENAME=honda96low.dat"
         STOP
      ENDIF
      
      DO 10 IZ=1,20
         READ(LUNHFLXL,'(A80)') DUMSTR
C         write(*,*) "Now reading ",DUMSTR
         DO 20 IE=1,72
            IF (IZ.le.11) THEN
               ICOS=IZ+9
            ELSE
               ICOS=21-IZ
            ENDIF
            READ(LUNHFLXL,*) HONEL_1D(IE),DUM(1),DUM(2),DUM(3),DUM(4),
     $                              HFLXLTMP(3),
     $                              HFLXLTMP(4),
     $                              HFLXLTMP(1),
     $                              HFLXLTMP(2),
     $                              DUM(5),DUM(6),DUM(7)
         HFLXL_1D(IE,ICOS,1,1)=HFLXLTMP(1)*HFLXLNRM_1D(1,1)
         HFLXL_1D(IE,ICOS,2,1)=HFLXLTMP(2)*HFLXLNRM_1D(2,1)
         HFLXL_1D(IE,ICOS,3,1)=HFLXLTMP(3)*HFLXLNRM_1D(3,1)
         HFLXL_1D(IE,ICOS,4,1)=HFLXLTMP(4)*HFLXLNRM_1D(4,1)

         HFLXL_1D(IE,ICOS,1,2)=HFLXLTMP(1)*HFLXLNRM_1D(1,2)
         HFLXL_1D(IE,ICOS,2,2)=HFLXLTMP(2)*HFLXLNRM_1D(2,2)
         HFLXL_1D(IE,ICOS,3,2)=HFLXLTMP(3)*HFLXLNRM_1D(3,2)
         HFLXL_1D(IE,ICOS,4,2)=HFLXLTMP(4)*HFLXLNRM_1D(4,2)
 20      CONTINUE
C      write(*,*) "IZ=",IZ," is OVER"
 10   CONTINUE

      IF ( IUSE_OPEN.EQ.0 ) THEN
         CALL SKCLOSEF(LUNHFLXL)
      ELSE
         CLOSE(LUNHFLXL)
      ENDIF

c
c     Rflist style open, kept for backwards compatibility
c
      IUSE_OPEN=0
      CALL SKOPENF(LUNHFLXM,1,'f',IHNDL)
      IF (IHNDL.LT.0) THEN
c
c        Check environment variable
c
         CALL GETENV("ATMPD_HONDAMID",FILENAME)
         IF (FILENAME(1:1) .EQ. ' ') THEN
c
c           Try for the default file with findconsts.
c
            CALL FINDCONSTS(FILENAME,'Honda Flux Mid','honda97mid.dat',2)
         END IF
         IF (FILENAME(1:1) .EQ. ' ') THEN
c
c           File still not found. Write error messages below
c           and exit just as before. (IHNDL still < 0)
c
            CONTINUE
         ELSE
            IHNDL = 0
            OPEN(LUNHFLXM,FILE=FILENAME,FORM='FORMATTED',STATUS='OLD')
            IUSE_OPEN=1
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
         write(*,*) "LUN=",LUNHFLXM,":FILENAME=honda97mid.dat"
         STOP
      ENDIF

      DO 40 IS=1,2
         DO 50 IZ=1,20
            DO 60 IA=1,12
               READ(LUNHFLXM,'(A80)') DUMSTR
C               write(*,*) "Now reading ",DUMSTR
               DO 70 IE=1,50
                  IF (IZ.le.11) THEN
                     ICOS=IZ+9
                  ELSE
                     ICOS=21-IZ
                  ENDIF
                  READ(LUNHFLXM,*) HONEM_1D(IE),
     $                                    HFLXM_1D(IE,ICOS,IA,3,IS),
     $                                    HFLXM_1D(IE,ICOS,IA,4,IS),
     $                                    HFLXM_1D(IE,ICOS,IA,1,IS),
     $                                    HFLXM_1D(IE,ICOS,IA,2,IS)
C                  write(*,'(E11.4,E11.4,E11.4,E11.4,E11.4)')
C     $                  HONEM(IE),
C     $                  HFLXM(IE,ICOS,IA,3,IS),
C     $                  HFLXM(IE,ICOS,IA,4,IS),
C     $                  HFLXM(IE,ICOS,IA,1,IS),
C     $                  HFLXM(IE,ICOS,IA,2,IS)
 70            CONTINUE
 60         CONTINUE
 50      CONTINUE
 40   CONTINUE

      IF ( IUSE_OPEN.EQ.0 ) THEN
         CALL SKCLOSEF(LUNHFLXM)
      ELSE
         CLOSE(LUNHFLXM)
      ENDIF

c
c     Rflist style open, kept for backwards compatibility
c
      IUSE_OPEN=0 
      CALL SKOPENF(LUNHFLXH,1,'f',IHNDL)
      IF (IHNDL.LT.0) THEN
c
c        Check environment variable
c
         CALL GETENV("ATMPD_HONDAHIGH",FILENAME)
         IF (FILENAME(1:1) .EQ. ' ') THEN
c
c           Try for the default file with findconsts.
c
            CALL FINDCONSTS(FILENAME,'Honda Flux High','honda96high.dat',2)
         END IF
         IF (FILENAME(1:1) .EQ. ' ') THEN
c
c           File still not found. Write error messages below
c           and exit just as before. (IHNDL still < 0)
c
            CONTINUE
         ELSE
            IHNDL = 0
            OPEN(LUNHFLXH,FILE=FILENAME,FORM='FORMATTED',STATUS='OLD')
            IUSE_OPEN=1 
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
         write(*,*) "LUN=",LUNHFLXH,":FILENAME=honda96high.dat"
         STOP
      ENDIF

      DO 80 IK=1,4
         READ(LUNHFLXH,'(A80)') DUMSTR
         DO 90 IE=1,31
            READ(LUNHFLXH,*)          HONEH_1D(IE)
     $           ,HFLXH_1D(IE,11,IK),HFLXH_1D(IE,12,IK)
     $           ,HFLXH_1D(IE,13,IK),HFLXH_1D(IE,14,IK)
     $           ,HFLXH_1D(IE,15,IK),HFLXH_1D(IE,16,IK)
     $           ,HFLXH_1D(IE,17,IK),HFLXH_1D(IE,18,IK)
     $           ,HFLXH_1D(IE,19,IK),HFLXH_1D(IE,20,IK)
            DO 100 IZ=1,10
               HFLXH_1D(IE,IZ+10,IK)
     $              =HFLXH_1D(IE,IZ+10,IK)/(HONEH_1D(IE)**3)
               HFLXH_1D(IE,11-IZ,IK)=HFLXH_1D(IE,IZ+10,IK)
 100        CONTINUE
 90      CONTINUE
 80   CONTINUE
      IF ( IUSE_OPEN.EQ.0 ) THEN
         CALL SKCLOSEF(LUNHFLXH)
      ELSE
         CLOSE(LUNHFLXH)
      ENDIF

      DO 200 IS=1,2
         DO 210 IK=1,4
            DO 220 IE=1,100
               HFLX_1D(IE,IK,IS)=0.
               IF (IE.le.35) THEN
                  HONE_1D(IE)=HONEL_1D(IE)
                  DO 230 IZ=1,20
                     HFLX_1D(IE,IK,IS)= HFLX_1D(IE,IK,IS)
     $                              +HFLXL_1D(IE,IZ,IK,IS)
 230              CONTINUE
                  HFLX_1D(IE,IK,IS)=HFLX_1D(IE,IK,IS)*0.1E0*2*3.141592
               ENDIF
               IF (IE.ge.36.and.IE.le.85) THEN
                  HONE_1D(IE)=HONEM_1D(IE-35)
                  DO 240 IZ=1,20
                     DO 250 IA=1,12
                        HFLX_1D(IE,IK,IS)= HFLX_1D(IE,IK,IS)
     $                                 +HFLXM_1D(IE-35,IZ,IA,IK,IS)
 250                 CONTINUE
 240              CONTINUE
                  HFLX_1D(IE,IK,IS)=HFLX_1D(IE,IK,IS)*0.1E0*(3.141592/6.)
               ENDIF
               IF (IE.ge.86) THEN
                  HONE_1D(IE)=HONEH_1D(IE-69)
                  DO 260 IZ=1,20
                     HFLX_1D(IE,IK,IS)=  HFLX_1D(IE,IK,IS)
     $                              +HFLXH_1D(IE-69,IZ,IK)
 260              CONTINUE
                  HFLX_1D(IE,IK,IS)=HFLX_1D(IE,IK,IS)*0.1E0*2*3.141592
               ENDIF
               HFLX_1D(IE,IK,IS)=HFLX_1D(IE,IK,IS)/(4.*3.141592)
 220        CONTINUE
 210     CONTINUE
 200  CONTINUE

      DO 300 IS=1,2
         DO 310 IK=1,4
            DO 320 IE=1,153
               HFLXA_1D(IE,IK,IS)=0.
               IF (IE.le.72) THEN
                  HONEA_1D(IE)=HONEL_1D(IE)
                  DO 330 IZ=1,20
                     HFLXA_1D(IE,IK,IS)= HFLXA_1D(IE,IK,IS)
     $                               +HFLXL_1D(IE,IZ,IK,IS)
 330              CONTINUE
                  HFLXA_1D(IE,IK,IS)=HFLXA_1D(IE,IK,IS)*0.1E0*2*3.141592
               ENDIF
               IF (IE.ge.73.and.IE.le.122) THEN
                  HONEA_1D(IE)=HONEM_1D(IE-72)
                  DO 340 IZ=1,20
                     DO 350 IA=1,12
                        HFLXA_1D(IE,IK,IS)= HFLXA_1D(IE,IK,IS)
     $                                 +HFLXM_1D(IE-72,IZ,IA,IK,IS)
 350                 CONTINUE
 340              CONTINUE
                  HFLXA_1D(IE,IK,IS)=HFLXA_1D(IE,IK,IS)*0.1E0*(3.141592/6.)
               ENDIF
               IF (IE.ge.123) THEN
                  HONEA_1D(IE)=HONEH_1D(IE-122)
                  DO 360 IZ=1,20
                     HFLXA_1D(IE,IK,IS)=  HFLXA_1D(IE,IK,IS)
     $                              +HFLXH_1D(IE-122,IZ,IK)
 360              CONTINUE
                  HFLXA_1D(IE,IK,IS)=HFLXA_1D(IE,IK,IS)*0.1E0*2*3.141592
               ENDIF
 320        CONTINUE
 310     CONTINUE
 300  CONTINUE


      END

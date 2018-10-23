************************************************************************
*     ---------------------------------
      SUBROUTINE NERDHON
*     ---------------------------------
C     
*     (Purpose)
C        READ FLUX CALCULATED BY M.Honda
C           
C     (Input) 
C        FLUX FOR LOW  ENERGY(  2MeV to 112MeV)/temporary KAMLOW
C        FLUX FOR MID  ENERGY(112MeV to  10TeV) /include high energy
C     (Creation Date and Author)
C       1996.??.?? ; Y.Hayato
C       2001.09.18 ; M.Ishitsuka 
C               replace flux table
C                  honda97mid.dat & honda96high.dat -> hkkm2.dat
C       2001.10.23 ; M.Ishitsuka
C               add solar activity to hkkm2.dat
C       2003.04.16 ; M.Ishitsuka
C               Update flux table to Honda2003 flux (hkkm03mt.dat)
C       2007.08.27 ; G.Mitsuka
C               Add Honda2006 flux (hkkm06mt.dat)
C       2007.11.16 ; G.Mitsuka
C               Add coeffecients to adjust Volkova flux to Honda flux
C       2007.12.05 ; G.Mitsuka
C               FC/PC's binning is merged with Upmu
C
************************************************************************

      IMPLICIT NONE
      
#include"hondaflx.h"
#include"necardap.h"

C      REAL*4 HFLXL(72,20,4,2),HFLXM(101,20,12,4,2),
C     $       HFLX(135,4,2),HFLXA(173,4,2)
C      REAL*4 HONE(135),HONEL(72),HONEM(101)
C      REAL*4 HONEA(72+101)
C      REAL*4 HFLXLNRM(4,2)
C      COMMON /HONFX/HONE ,HFLX,
C     $              HONEA,HFLXA,
C     $              HONEL,HFLXL,
C     $              HONEM,HFLXM,
C     $              HFLXLNRM

      INTEGER*4 LUNHFLXL
      PARAMETER (LUNHFLXL=77)
      INTEGER*4 LUNHFLXM
c      PARAMETER (LUNHFLXM=78)

      INTEGER*4 IHNDL

      INTEGER*4 IK,IE,IZ,IA,ICOS,IS
      CHARACTER*256  FILENAME
      CHARACTER*80   DUMSTR
      REAL*4    DUM(10),HFLXLTMP(4)

      INTEGER*4 ILOADED
      SAVE      ILOADED
      DATA      ILOADED/0/

      INTEGER*4 IUSE_OPEN

***** above .100E+05 GeV is volkova shape *****
      INTEGER*4 I,J,ICOS_UM
      REAL*4 ENE_UM(8),FNUM1(10,8),FNUM2(10,8)
      
      DATA ENE_UM/.316D+04,.100E+05,.300E+05,.100E+06,
     $            .300E+06,.100E+07,.300E+07,.100E+08/
*      DATA COS_UM/ .05D+00,.15D+00,.25D+00,.35D+00,.45D+00,
*     &             .55D+00,.65D+00,.75D+00,.85D+00,.95D+00/
* nu_mu * E^3 [GeV^2/(m^2 s sr)]
      DATA ((FNUM1(I,J),I=1,10),J=1,8)/
     +413.8,289.0,216.2,173.8,150.5,138.7,129.2,120.4,112.5,105.2,
     +0.3002D+03, 0.1987D+03, 0.1382D+03, 0.1046D+03, 0.8675D+02,
     +0.7715D+02, 0.6979D+02, 0.6350D+02, 0.5814D+02, 0.5340D+02,
     &0.1816D+03, 0.1151D+03, 0.7684D+02, 0.5652D+02, 0.4583D+02,
     &0.4016D+02, 0.3594D+02, 0.3240D+02, 0.2944D+02, 0.2685D+02,
     &0.8730D+02, 0.5288D+02, 0.3421D+02, 0.2475D+02, 0.1989D+02,
     &0.1737D+02, 0.1553D+02, 0.1400D+02, 0.1271D+02, 0.1159D+02,
     &0.4445D+02, 0.2662D+02, 0.1690D+02, 0.1214D+02, 0.9766D+01,
     &0.8452D+01, 0.7490D+01, 0.6725D+01, 0.6111D+01, 0.5588D+01,
     &0.1578D+02, 0.9449D+01, 0.6002D+01, 0.4351D+01, 0.3525D+01,
     &0.3023D+01, 0.2654D+01, 0.2377D+01, 0.2165D+01, 0.1992D+01,
     &0.5154D+01, 0.3082D+01, 0.1955D+01, 0.1417D+01, 0.1146D+01,
     &0.9835D+00, 0.8633D+00, 0.7727D+00, 0.7033D+00, 0.6460D+00,
     &0.1486D+01, 0.8884D+00, 0.5658D+00, 0.4084D+00, 0.3282D+00,
     &0.2840D+00, 0.2517D+00, 0.2260D+00, 0.2054D+00, 0.1878D+00/
* nu_mubar * E^3 [GeV^2/(m^2 s sr)]
      DATA ((FNUM2(I,J),I=1,10),J=1,8)/
     +271.0,218.3,177.2,145.5,121.0,101.5, 85.4, 72.0, 60.7, 51.3,
     +0.1819D+03, 0.1389D+03, 0.1049D+03, 0.8103D+02, 0.6455D+02,
     +0.5225D+02, 0.4270D+02, 0.3514D+02, 0.2903D+02, 0.2410D+02,
     &0.1028D+03, 0.7514D+02, 0.5442D+02, 0.4089D+02, 0.3184D+02,
     &0.2539D+02, 0.2053D+02, 0.1674D+02, 0.1373D+02, 0.1131D+02,
     &0.4606D+02, 0.3217D+02, 0.2259D+02, 0.1669D+02, 0.1288D+02,
     &0.1024D+02, 0.8269D+01, 0.6742D+01, 0.5525D+01, 0.4553D+01,
     &0.2208D+02, 0.1525D+02, 0.1051D+02, 0.7711D+01, 0.5956D+01,
     &0.4692D+01, 0.3756D+01, 0.3051D+01, 0.2501D+01, 0.2067D+01,
     &0.7371D+01, 0.5090D+01, 0.3508D+01, 0.2597D+01, 0.2021D+01,
     &0.1578D+01, 0.1251D+01, 0.1013D+01, 0.8331D+00, 0.6925D+00,
     &0.2282D+01, 0.1574D+01, 0.1083D+01, 0.8017D+00, 0.6231D+00,
     &0.4866D+00, 0.3858D+00, 0.3124D+00, 0.2565D+00, 0.2130D+00,
     &0.6226D+00, 0.4292D+00, 0.2966D+00, 0.2187D+00, 0.1688D+00,
     &0.1329D+00, 0.1064D+00, 0.8647D-01, 0.7090D-01, 0.5859D-01/

*     Scaling factor between Honda flux and Volkova flux
*     Coeffecients are dereved by fitting a flux around 10TeV
*     (G.Mitsuka, 2007 Nov. 16th)
      REAL VOL_COEFF(4)

      if (NEFLXFLG.eq.0) then   ! Honda06
         VOL_COEFF(1) = 0.03039
         VOL_COEFF(2) = 0.03375
         VOL_COEFF(3) = 0.96184
         VOL_COEFF(4) = 0.84404
      else                      ! Honda03 etc.
         VOL_COEFF(1) = 1.6E-2
         VOL_COEFF(2) = 2.1E-2
         VOL_COEFF(3) = 5.4E-1
         VOL_COEFF(4) = 5.7E-1
      endif

C---- Set LUNHFLXM by NEFLXFLG
      if(NEFLXFLG.eq.0) then
         LUNHFLXM = 78
      elseif(NEFLXFLG.eq.6) then
         LUNHFLXM = 79
      elseif(NEFLXFLG.lt.0.or.NEFLXFLG.gt.7) then
         print*, "Invalid NEFLXFLG value : nerdhon.F"
         stop
      endif

C---- CHECK CALLED TWICE?
      IF (ILOADED.ne.0) RETURN
      if(NEFLXFLG.eq.0) write(*,*) "Now LOADING HONDA 2006 FLUX IN NEUTFLUX"
      if(NEFLXFLG.eq.6) write(*,*) "Now LOADING HONDA 2003 FLUX IN NEUTFLUX"
      ILOADED = 1

C---  Normalize Factor for L.E.Table
C     Temporary(To adjust old-table and new table)
C     Factor was decided to match total flax

      HFLXLNRM(1,1)=1.00
      HFLXLNRM(2,1)=1.05
      HFLXLNRM(3,1)=1.00
      HFLXLNRM(4,1)=1.00

      HFLXLNRM(1,2)=0.89
      HFLXLNRM(2,2)=0.95
      HFLXLNRM(3,2)=0.89
      HFLXLNRM(4,2)=0.89
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
            READ(LUNHFLXL,*) HONEL(IE),DUM(1),DUM(2),DUM(3),DUM(4),
     $                              HFLXLTMP(3),
     $                              HFLXLTMP(4),
     $                              HFLXLTMP(1),
     $                              HFLXLTMP(2),
     $                              DUM(5),DUM(6),DUM(7)
         HFLXL(IE,ICOS,1,1)=HFLXLTMP(1)*HFLXLNRM(1,1)
         HFLXL(IE,ICOS,2,1)=HFLXLTMP(2)*HFLXLNRM(2,1)
         HFLXL(IE,ICOS,3,1)=HFLXLTMP(3)*HFLXLNRM(3,1)
         HFLXL(IE,ICOS,4,1)=HFLXLTMP(4)*HFLXLNRM(4,1)

         HFLXL(IE,ICOS,1,2)=HFLXLTMP(1)*HFLXLNRM(1,2)
         HFLXL(IE,ICOS,2,2)=HFLXLTMP(2)*HFLXLNRM(2,2)
         HFLXL(IE,ICOS,3,2)=HFLXLTMP(3)*HFLXLNRM(3,2)
         HFLXL(IE,ICOS,4,2)=HFLXLTMP(4)*HFLXLNRM(4,2)

 20      CONTINUE
C      write(*,*) "IZ=",IZ," is OVER"
 10   CONTINUE

      IF ( IUSE_OPEN.eq.0 ) THEN
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
         if(NEFLXFLG.eq.0) then
            CALL FINDCONSTS(FILENAME,'Honda Flux Mid & High','hkkm06mt.dat',2)
         elseif(NEFLXFLG.eq.6) then
            CALL FINDCONSTS(FILENAME,'Honda Flux Mid & High','hkkm03mt.dat',2)
         endif
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
         if(NEFLXFLG.eq.0) write(*,*) "LUN=",LUNHFLXM,":FILENAME=hkkm06mt.dat"
         if(NEFLXFLG.eq.6) write(*,*) "LUN=",LUNHFLXM,":FILENAME=hkkm03mt.dat"
         STOP
      ENDIF

      DO 40 IS=1,2
         DO 50 IZ=1,20
            DO 60 IA=1,12
               READ(LUNHFLXM,'(A80)') DUMSTR
               READ(LUNHFLXM,'(A80)') DUMSTR
               DO 70 IE=1,101
                  ICOS=21-IZ
                  READ(LUNHFLXM,*) HONEM(IE),
     $                 HFLXM(IE,ICOS,IA,3,IS),
     $                 HFLXM(IE,ICOS,IA,4,IS),
     $                 HFLXM(IE,ICOS,IA,1,IS),
     $                 HFLXM(IE,ICOS,IA,2,IS)
 70            CONTINUE
 60         CONTINUE
 50      CONTINUE
 40   CONTINUE

      IF ( IUSE_OPEN.eq.0 ) THEN
         CALL SKCLOSEF(LUNHFLXM)
      ELSE
         CLOSE(LUNHFLXM)
      ENDIF

      DO 140 IS=1,2   !--- solact
         DO 150 ICOS=1,20   !--- zenith angle
            DO 160 IA=1,12  !--- azimith angle
#ifdef UPMU_RCK
               DO 170 IE=100,107  !--- energy
#else
               DO 170 IE=102,107  !--- energy
#endif
                  IF (ICOS.le.10) THEN
                     ICOS_UM=11-ICOS
                  ELSEIF(ICOS.GE.11) THEN
                     ICOS_UM=21-ICOS
                  ELSE
                     STOP
                  ENDIF
                  HONEM(IE)=ENE_UM(IE-99)
                  HFLXM(IE,ICOS,IA,3,IS)=
     $                 FNUM1(ICOS_UM,IE-99)/(HONEM(IE)**3)*VOL_COEFF(3)
                  HFLXM(IE,ICOS,IA,4,IS)= 
     $                 FNUM2(ICOS_UM,IE-99)/(HONEM(IE)**3)*VOL_COEFF(4)
                  HFLXM(IE,ICOS,IA,1,IS)=
     $                 FNUM1(ICOS_UM,IE-99)/(HONEM(IE)**3)*VOL_COEFF(1)
                  HFLXM(IE,ICOS,IA,2,IS)=
     $                 FNUM2(ICOS_UM,IE-99)/(HONEM(IE)**3)*VOL_COEFF(2)
 170           CONTINUE
 160        CONTINUE
 150     CONTINUE
 140  CONTINUE


      DO 200 IS=1,2
         DO 210 IK=1,4
            DO 220 IE=1,141
               HFLX(IE,IK,IS)=0.
               IF (IE.le.34) THEN
                  HONE(IE)=HONEL(IE)
                  DO 230 IZ=1,20
                     HFLX(IE,IK,IS)= HFLX(IE,IK,IS)
     $                              +HFLXL(IE,IZ,IK,IS)
 230              CONTINUE
                  HFLX(IE,IK,IS)=HFLX(IE,IK,IS)*0.1E0*2*3.141592
               ENDIF
               IF (IE.ge.35) THEN
                  HONE(IE)=HONEM(IE-34)
                  DO 240 IZ=1,20
                     DO 250 IA=1,12
                        HFLX(IE,IK,IS)= HFLX(IE,IK,IS)
     $                                 +HFLXM(IE-34,IZ,IA,IK,IS)
 250                 CONTINUE
 240              CONTINUE
                  HFLX(IE,IK,IS)=HFLX(IE,IK,IS)*0.1E0*(3.141592/6.)
               ENDIF
               HFLX(IE,IK,IS)=HFLX(IE,IK,IS)/(4.*3.141592)
 220        CONTINUE
 210     CONTINUE
 200  CONTINUE

      DO 300 IS=1,2
         DO 310 IK=1,4
            DO 320 IE=1,173
               HFLXA(IE,IK,IS)=0.
               IF (IE.le.72) THEN
                  HONEA(IE)=HONEL(IE)
                  DO 330 IZ=1,20
                     HFLXA(IE,IK,IS)= HFLXA(IE,IK,IS)
     $                               +HFLXL(IE,IZ,IK,IS)
 330              CONTINUE
                  HFLXA(IE,IK,IS)=HFLXA(IE,IK,IS)*0.1E0*2*3.141592
               ENDIF
               IF (IE.ge.73) THEN
                  HONEA(IE)=HONEM(IE-72)
                  DO 340 IZ=1,20
                     DO 350 IA=1,12
                        HFLXA(IE,IK,IS)= HFLXA(IE,IK,IS)
     $                                 +HFLXM(IE-72,IZ,IA,IK,IS)
 350                 CONTINUE
 340              CONTINUE
                  HFLXA(IE,IK,IS)=HFLXA(IE,IK,IS)*0.1E0*(3.141592/6.)
               ENDIF
 320        CONTINUE
 310     CONTINUE
 300  CONTINUE

      END

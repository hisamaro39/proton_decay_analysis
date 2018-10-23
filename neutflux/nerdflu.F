************************************************************************
*     ---------------------------------
      SUBROUTINE NERDFLU
*     ---------------------------------
C     
*     (Purpose)
C        READ FLUX CALCULATED BY FLUKA
C           
C     (Input) 
C        differential FLUX 
C
C     (Creation Date and Author)
C        2000.10.13 ; A. Kibayashi
C
************************************************************************

      IMPLICIT NONE
      
#include"flukaflx.h"

      INTEGER*4 LUNFFLX 
      PARAMETER (LUNFFLX=66)
      INTEGER*4 LUNFFLXMX
      PARAMETER (LUNFFLXMX=67)   
      INTEGER*4 IFLKL

      INTEGER*4 IK,IE,IZ,ICOS,IS
      CHARACTER*80  DUMSTR
      CHARACTER*256 FILENAME
c      REAL*4    DUM(1),FFLXTMP(4)
      REAL*4    DUM(1)

      INTEGER*4 ILOADED
      SAVE      ILOADED
      DATA      ILOADED/0/

      
C---- CHECK CALLED TWICE?
      IF (ILOADED.ne.0) RETURN
      write(*,*) "Now LOADING FLUKA FLUX"
      ILOADED = 1


*-ak fluka flux at solar min.
c
c     Rflist style open, kept for backwards compatibility
c
      CALL SKOPENF(LUNFFLX,1,'f',IFLKL)
      IF (IFLKL.LT.0) THEN
c
c        Check environment variable
c
         CALL GETENV("ATMPD_FLUKA",FILENAME)
         IF (FILENAME(1:1) .EQ. ' ') THEN
c
c           Try for the default file with findconsts.
c
            CALL FINDCONSTS(FILENAME,'Fluka Flux','fit_fluka02.dat',2)
         END IF
         IF (FILENAME(1:1) .EQ. ' ') THEN
c
c           File still not found. Write error messages below
c           and exit just as before. (IFLKL still < 0)
c
            CONTINUE
         ELSE
            IFLKL = 0
            OPEN(LUNFFLX,FILE=FILENAME,FORM='FORMATTED',STATUS='OLD')
         END IF
      END IF
c
c     Old RFLIST behavior- kept for backwards compatibility...
c
      IF (IFLKL.LT.0) THEN
         IF (IFLKL.eq.-1001) THEN
            write(*,*) "******    ERROR    ******"
            write(*,*) "You did not specify the flux file!"
         ENDIF
         write(*,*) "******    ERROR    ******"
         write(*,*) "Could not open the FLUX file!"
         write(*,*) "LUN=",LUNFFLX,":FILENAME=fit_fluka02.dat"
         STOP
      ENDIF
      
      DO 10 IK=1,4              ! neutrino flavour
         READ(LUNFFLX,'(A80)') DUMSTR
c         write(*,*) "Now reading ",DUMSTR
         DO 20 IZ=1,40          ! cos zenith
            ICOS=IZ
            do 30 IE=1,61       ! neutrino energy
               READ(LUNFFLX,900) FONE(IE),DUM(1),FFLX(IE,ICOS,IK,1)

 30         continue
 20      CONTINUE
c      write(*,*) "IZ=",IZ," is OVER"
 10   CONTINUE

      CALL SKCLOSEF(LUNFFLX)

*-ak fluka flux at solar max.
c
c     Rflist style open, kept for backwards compatibility
c
      CALL SKOPENF(LUNFFLXMX,1,'f',IFLKL)
      IF (IFLKL.LT.0) THEN
c
c        Check environment variable
c
         CALL GETENV("ATMPD_FLUKAMX",FILENAME)
         IF (FILENAME(1:1) .EQ. ' ') THEN
c
c           Try for the default file with findconsts.
c
            CALL FINDCONSTS(FILENAME,'Fluka Flux','fit_fluka02_mx.dat',2)
         END IF
         IF (FILENAME(1:1) .EQ. ' ') THEN
c
c           File still not found. Write error messages below
c           and exit just as before. (IFLKL still < 0)
c
            CONTINUE
         ELSE
            IFLKL = 0
            OPEN(LUNFFLXMX,FILE=FILENAME,FORM='FORMATTED',STATUS='OLD')
         END IF
      END IF
c
c     Old RFLIST behavior- kept for backwards compatibility...
c
      IF (IFLKL.LT.0) THEN
         IF (IFLKL.eq.-1001) THEN
            write(*,*) "******    ERROR    ******"
            write(*,*) "You did not specify the flux file!"
         ENDIF
         write(*,*) "******    ERROR    ******"
         write(*,*) "Could not open the FLUX file!"
         write(*,*) "LUN=",LUNFFLXMX,":FILENAME=fit_fluka02_mx.dat"
         STOP
      ENDIF
      
      DO 40 IK=1,4              ! neutrino flavour
         READ(LUNFFLXMX,'(A80)') DUMSTR
c         write(*,*) "Now reading ",DUMSTR
         DO 50 IZ=1,40          ! cos zenith
            ICOS=IZ
            do 60 IE=1,61       ! neutrino energy
               READ(LUNFFLXMX,900) FONE(IE),DUM(1),FFLX(IE,ICOS,IK,2)

 60         continue
 50      CONTINUE
c      write(*,*) "IZ=",IZ," is OVER"
 40   CONTINUE



*--- integration ---*
      DO 200 IS=1,2             ! high and low for solar act.
         DO 210 IK=1,4          ! neutrino flavour
            DO 220 IE=1,61      ! neutrino energy
               FFLXA(IE,IK,IS)=0.
                  FONEA(IE)=FONE(IE)
                  DO 230 IZ=1,40 ! sum all cosz directions
                     FFLXA(IE,IK,IS)= FFLXA(IE,IK,IS)
     $                              +FFLX(IE,IZ,IK,IS)
 230              CONTINUE
                  FFLXA(IE,IK,IS)=FFLXA(IE,IK,IS)*0.1E0*3.141592
 220        CONTINUE
 210     CONTINUE
 200  CONTINUE

 900  FORMAT (1x,f10.3,3x,f10.3,3x,f20.8)

      END


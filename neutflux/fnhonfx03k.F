      SUBROUTINE FNHONFX03K(ENUE,ORGDIR,IPAR,FLX)

c Return differential Honda flux GeV-1/m^2/s/sr
c
c Input:
c     ENUE    Energy of neutrino (GeV)
c     ORGDIR  Direction of neutrino origin (Super-K convention)
c             orgdir(3)<0 : downward direction
c             orgdir(3)>0 : upward direction
c     IPAR    Flavor of neutrino(12=nue,-12=nuebar,14=numu,-14=numubar)
c
c     Mar 2004    Masaki Ishitsuka 
c                 based on fnhonfx.F (Y.Hayato)
c
c     format of the table : hkkm03mt.dat
c ============================================================================
c  unit    /GeV^1/m^2/s/sr
c ============================================================================
c  comment (cos=*)
c  number  (energy nue sig_nue nuebar sig_nuebar num sig_num numbar sig_numbar)
c ============================================================================
      IMPLICIT NONE

      real flx(2)

      INTEGER*4 LUNHFLX
      PARAMETER (LUNHFLX=65)

      INTEGER*4 IHNDL

      REAL ENERGY(118)
      REAL BFNU(118,10,4)
      REAL BFNUK(118,10,4)
c ----------------------------------------------------------------------
      REAL FINT
      EXTERNAL FINT
      REAL FNHONFX
      EXTERNAL FNHONFX
      CHARACTER*132 DUMMY_STRING, FILENAME
      REAL ENUE, ORGDIR(3), DIR(3)
      INTEGER IPAR 
      INTEGER J
      REAL BIN
      INTEGER ICOSZEN, ILOGE
      LOGICAL INITIALIZED
      DATA INITIALIZED/.FALSE./
      INTEGER*4  IEL,IEH
      INTEGER*4  ITL,ITH
      INTEGER*4  IPL,IPH

      REAL*4     RTL,RTH
      REAL*4     RPL,RPH

      REAL*4     FXLELTPL,FXLELTPH,FXHELTPL,FXHELTPH
      REAL*4     FXLEHTPL,FXLEHTPH,FXHEHTPL,FXHEHTPH
      REAL*4     FXLELTPLK,FXLELTPHK,FXHELTPLK,FXHELTPHK
      REAL*4     FXLEHTPLK,FXLEHTPHK,FXHEHTPLK,FXHEHTPHK
      
      REAL*4     FLXLELT,FLXLEHT
      REAL*4     FLXHELT,FLXHEHT
      REAL*4     FLXLE,FLXHE
      REAL*4     FLXLELTK,FLXLEHTK
      REAL*4     FLXHELTK,FLXHEHTK
      REAL*4     FLXLEK,FLXHEK
      REAL*4     PHI,APHI

      REAL*4     A,B
      REAL*4     AK,BK

      INTEGER*4  IPKIND
      real dum1(2), dum2(2)
      
c ======================================================================
c
c     Initialize upon first call by loading tables:
c     =============================================
      IF (.NOT.INITIALIZED) THEN
         WRITE(*,'(A)') 
     &        " Now LOADING RVW HONDA-2003 FLUX for K/Pi ratio  "
c
c        Load solar min file:
c        --------------------
         CALL SKOPENF(LUNHFLX,1,'f',IHNDL)
         IF (IHNDL.LT.0) THEN
            print *, 'Could not open hkkm03-kaon flux file'
            stop
         ENDIF
c         CALL FINDCONSTS(FILENAME,'Honda Flux 2003',
c     &        'hkkm03.dat',0)
c         OPEN(65,FILE=FILENAME,FORM='FORMATTED',STATUS='OLD')
c     ------------------zenith-----------------------------------
         DO ICOSZEN = 1, 10
c     ------------------energy-----------------------------------
            DO ILOGE = 1, 118
               READ(LUNHFLX,*) ENERGY(ILOGE),
     &              BFNU(ILOGE,11-ICOSZEN,3),
     &              BFNU(ILOGE,11-ICOSZEN,4),
     &              BFNU(ILOGE,11-ICOSZEN,1),
     &              BFNU(ILOGE,11-ICOSZEN,2),
     &              dum1,
     &              BFNUK(ILOGE,11-ICOSZEN,3),
     &              BFNUK(ILOGE,11-ICOSZEN,4),
     &              BFNUK(ILOGE,11-ICOSZEN,1),
     &              BFNUK(ILOGE,11-ICOSZEN,2),
     &              dum2
            END DO
         END DO
         CALL SKCLOSEF(LUNHFLX)
         INITIALIZED = .TRUE.
      END IF                    ! End of first-time initialization

C---CHECK KIND
      IF (IPAR.eq.12)  IPKIND=1
      IF (IPAR.eq.-12) IPKIND=2
      IF (IPAR.eq.14)  IPKIND=3
      IF (IPAR.eq.-14) IPKIND=4

c
c     Geometry conversion, made to agree with neut/fnhonfx.F
c     ================================================
      DIR(3) = -ORGDIR(3)
      DIR(1)=-(COS(3.141593*2435./60./180.)*ORGDIR(1)
     &       +SIN(3.141593*2435./60./180.)*ORGDIR(2))
      DIR(2)=-(-SIN(3.141593*2435./60./180.)*ORGDIR(1)
     &       +COS(3.141593*2435./60./180.)*ORGDIR(2))
c
c     At high energy or low energy, use honda-2001 routine:
c     ================================
      IF ((ENUE.GT.ENERGY(118)).or.(ENUE.LT.ENERGY(1))) THEN
         FLX(1) = FNHONFX(ENUE,ORGDIR,IPAR)
         FLX(2) = 0.
      ELSE IF ((ENUE.LE.ENERGY(118)) .AND. (ENUE.GE.ENERGY(1))) THEN
         DO 5 IEH=1,118
            IF (ENUE.lt.ENERGY(IEH)) GOTO 7
 5       CONTINUE
 7       IEL=IEH-1

         IF (DIR(3).le.-0.95) THEN
            ITL=1
            ITH=1
            GOTO 10
         ENDIF
         IF (DIR(3).gt.-0.05) THEN
            ITL=10
            ITH=10
            GOTO 10
         ENDIF
         ITL=INT(DIR(3)*10.+10.5)
         ITH=ITL+1
         IF (ITH.eq.11) THEN 
            ITH=10
            GOTO 10
         ENDIF
         
 10      CONTINUE
         
C---CONSIDER EACH FLUX for LOW ENERGY AND HIGH ENERGY

 15      CONTINUE

         RTL     = REAL(ITL)*0.1-1.05
         RTH     = REAL(ITH)*0.1-1.05
         
         FLXLELT = BFNU(IEL,ITL,IPKIND)
         FLXLEHT = BFNU(IEL,ITH,IPKIND)
         FLXLELTK = BFNUK(IEL,ITL,IPKIND)
         FLXLEHTK = BFNUK(IEL,ITH,IPKIND)
         
         IF (ITL.eq.ITH) THEN
            FLXLE=FLXLELT
            FLXLEK=FLXLELTK
            GOTO 20
         ENDIF
         
         FLXLE   = (FLXLEHT-FLXLELT)/(RTH-RTL)*(DIR(3)-RTL)+FLXLELT
         FLXLEK   = (FLXLEHTK-FLXLELTK)/(RTH-RTL)*(DIR(3)-RTL)+FLXLELTK
         
 20      FLXHELT  = BFNU(IEH,ITL,IPKIND)
         FLXHEHT  = BFNU(IEH,ITH,IPKIND)
         FLXHELTK  = BFNUK(IEH,ITL,IPKIND)
         FLXHEHTK  = BFNUK(IEH,ITH,IPKIND)
         
         IF (ITL.eq.ITH) THEN
            FLXHE=FLXHELT
            FLXHEK=FLXHELTK
            GOTO 30
         ENDIF
         
         FLXHE   = (FLXHEHT-FLXHELT)/(RTH-RTL)*(DIR(3)-RTL)+FLXHELT
         FLXHEK   = (FLXHEHTK-FLXHELTK)/(RTH-RTL)*(DIR(3)-RTL)+FLXHELTK
         
 30      B       = ALOG(FLXLE/FLXHE)/ALOG(ENERGY(IEL)/ENERGY(IEH))
         A       = FLXLE/ENERGY(IEL)**B
         FLX(1) = A*ENUE**B
         BK       = ALOG(FLXLEK/FLXHEK)/ALOG(ENERGY(IEL)/ENERGY(IEH))
         AK       = FLXLEK/ENERGY(IEL)**BK
         FLX(2) = AK*ENUE**BK

      END IF

      RETURN
      END

************************************************************************
*     ----------------------------------------------
      SUBROUTINE NEEVRATE(NSOL, SACT, EVRATE, RNEUT)
*     ----------------------------------------------
*
*     (Purpose)
*       Atmospehric neutrino event rate / kton yr
*
*     (Input)
*       NSOL   : FLAG OF SOLAR ACTIVETY
*                0 : SACT = 0.5         
*       SACT   : SOLAR ACTIVETY
*
*     (Output)
*       EVRATE   : EVENT RATE / KTON YR
*       RNEUT(4) : RATIO OF EACH NEUTRINO SPECIES 
*                  I=1 RATIO OF NEU(E)
*                  I=2 RATIO OF NEU(E)+NEU(E)BAR
*                  I=3 RATIO OF NEU(E)+NEU(E)BAR+NEU(MU)
*                  I=4 RATIO OF NEU(E)+NEU(E)BAR+NEU(MU)+NEU(MU)BAR=1.
*
*     (Creation Date and Author)
C       1988-11-26 ; T. KAJITA
*       1989-07-18 ; T. KAJITA
C                       (FINAL VERSION)
*       1994.10.31 ; K. KANEYUKI  FOR S.K.
*                    EVRATE -> NEEVRATE
C       1994-09-03 ; T. KAJITA
C                      (DATA ARE AUTOMATICALLY UPDATED WHENEVER
C                       THE FLUX AND THE CROSS SECTIONS ARE UPDATED)
C       1995-05-18 ; Y. HAYATO
C                      ORIGINAL LATEST VERSION FOR S.K.
C       2004-02-04 ; C.Saji
C                      ADD TAU MODE
C       2007-11-10 ; T.Tanaka
*                      ADD UPMU MODE
C       2007-12-05 ; G.Mitsuka
C                      FC/PC's binning is merged with Upmu
************************************************************************
#include <nesolact.h>
#include "necardap.h"
#include "necard.h"       


#ifdef UPMU_RCK
#include "upmu.h"
#endif

C      REAL RNEUT(4),RIPMAX(4),RIPMIN(4),RIP(4)
      REAL RNEUT(4),RIP(4)
      REAL FLUXVL

C            IF(SOLACT=1.) -- SOLAR MAX.
C            IF(SOLACT=0.) -- SOLAR MIN.
C         (  IF(SOLACT=0.5) -- AVERAGE OF SOLAR MIN/MAX )
C         (  IF(NSOL=0) -- SOLACT=0.5 UNCONDITIONALLY )
C
C
#ifndef UPMU_RCK
      DATA PI/3.141593/
#endif

      DATA PNUM/6.03E32/
C              TOTAL NUMBER OF NUCLEON IN 1000TON OF WATER
      DATA TIME/3.15576E7/
C              TOTAL SECOND IN A YEAR
C
C
      ISOL = NSOL
      SOLACT = SACT

      IF(ISOL.EQ.0) SOLACT=0.5
      IF(SOLACT.LT.0.) SOLACT=0.
      IF(SOLACT.GT.1.) SOLACT=1.
C
C  CALCULATION FROM #LOWHIGH
C
      DO 100 I=1,4
         RIP(I)=0.0
         IF(I.EQ.1) IPAR=12
         IF(I.EQ.2) IPAR=-12
         IF(ITAUFLG.EQ.1) THEN
            IF(I.EQ.3)IPAR=16
            IF(I.EQ.4)IPAR=-16
         ELSE
            IF(I.EQ.3)IPAR=14
            IF(I.EQ.4)IPAR=-14
         ENDIF
C     


#ifdef UPMU_RCK
*--- Minimum neutrino energy (GeV)
      emax=UPMU_ETH
*---->  Energy mx up to 100TeV(C.Saji 00/09/29)
      DO 101 III=1,5
         RIPTMP=0.D0
         NINT=100
         emin=emax
         if(iii.eq.1) emax=10.
         if(iii.eq.2) emax=100.
         if(iii.eq.3) emax=1000.
*------->   energy up to 100TeV(C.Saji 00/09/29)
         if(iii.eq.4) emax=10000.
         if(iii.eq.5) emax=90000.
#else

         DO 101 III=1,7
            RIPTMP=0.D0
            NINT=100
            if (III.eq.1) then     ! 10MeV  to  200MeV
               EMIN=0.01
               EMAX=0.2
            elseif (III.eq.2) then ! 200MeV to  5GeV
               EMIN=0.2
               EMAX=5.
            elseif (III.eq.3) then ! 5GeV   to  10GeV
               EMIN=5.
               EMAX=10.
            elseif (III.eq.4) then ! 10GeV  to  100GeV
               EMIN=10.
               EMAX=100.
            elseif (III.eq.5) then ! 100GeV to  1TeV
               EMIN=100.
               EMAX=1000.
            elseif (III.eq.6) then ! 1TeV   to  10TeV
               EMIN=1000.
               EMAX=10000.
            elseif (III.eq.7) then ! 10TeV  to  90TeV
               EMIN=10000.
               EMAX=90000.
            endif

C     300Steps only between 5GeV and 1TeV for the compatibility with the old binning
            IF(III.GE.3.AND.III.LE.5) NINT=300

#endif
            EBIN=(EMAX-EMIN)/FLOAT(NINT)
C     
C     ++ NEW INTEGRAL
C     
               DO 200 J=1,NINT
               E=EMIN+EBIN*(FLOAT(J)-0.5)
               CALL NEFLUX(IPAR,E,FLUXVL)
C               write(*,*) IPAR,":",E,":",FLUXVL
C               RIPTMP=RIPTMP+FLUXVL*FNTOTPAU(IPAR,E)*1e-38
               RIPTMP=RIPTMP+FLUXVL*FNTOTPAU(IPAR,E)

C       write(*,*) "FLUXVL=",FLUXVL," FNTOTPAU(",IPAR,",",E,")=",
C     $                     FNTOTPAU(IPAR,E)

 200        CONTINUE
            RIPTMP=RIPTMP*1e-30
            RIPTMP=RIPTMP*1e-8
C            RIP(I)=RIP(I)+RIPTMP*EBIN*1.E-4*PNUM*4.*PI*TIME
            RIP(I)=RIP(I)+RIPTMP*PNUM*TIME*EBIN*1.E-4*4.*PI
C           write(*,*) "RIP(",i,")=",RIP(I)
 101     CONTINUE
C     
 100  CONTINUE
C     
C     
      EVRATE=0.
      DO 30 I=1,4
         EVRATE=EVRATE+RIP(I)
 30   CONTINUE
C     
      RNEUT(1)=RIP(1)/EVRATE
      RNEUT(2)=(RIP(1)+RIP(2))/EVRATE
      RNEUT(3)=(RIP(1)+RIP(2)+RIP(3))/EVRATE
      RNEUT(4)=(RIP(1)+RIP(2)+RIP(3)+RIP(4))/EVRATE
      IF(RNEUT(4).LT.0.999.OR.RNEUT(4).GT.1.001) THEN
         WRITE(6,600) RNEUT(4)
 600     FORMAT(' ',' ************* ERROR RNEUT(4)= ',F7.4,' *********')
      ELSE
         RNEUT(4)=1.
      END IF
C     
      RETURN
      END

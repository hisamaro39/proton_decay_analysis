************************************************************************
*     ---------------------------------
      SUBROUTINE NEMAKE(LUNO,YEAR)
*     ---------------------------------
C     
*     (Purpose)
C        Generate atmospheric neutrino events.
C      
C     (Input)
C       LUNO   : Logical Unit Number for output
C       YEAR   : YEAR
C
*     (Output)
C      
C     (Creation Date and Author)
C       1996.??.?? ; Y.Hayato
C       2007.11.10 ; T.Tanaka
C                    ADD UPMU MODE                         
C       2007.12.05 ; G.Mitsuka
C                    Simulated volume of FC/PC is expanded 
C                     to the entire water tank
************************************************************************
      IMPLICIT NONE

#include "dsbnkf.h"
#include "necardap.h"
#include "necardvc.h"
#include "nesolact.h"
#include "necard.h"

#ifdef UPMU_RCK
#include "upmu.h"
#include "geotnk.h"
      REAL*4 CYLVKTON           !--cylinder volkton
      REAL*4 UPMU_EVFACT        !--event volume factor
      INTEGER*4 STORE
      INTEGER*4 UFLAG
#endif

      INTEGER*4 LUNO
      REAL*4    YEAR,DWALL

      REAL*4    VOLKTON,EVRATE

      INTEGER*4 IPAR,NUMEVT,NUMREM,I,IERR
      
      INTEGER*4 ISOLR
      REAL*4    SOLACTV

      INTEGER*4 LUX,ISEED1,ISEED2
      INTEGER*4 K1,K2,K3,K4
      INTEGER*4 ISTAU1,ISTAU2,ISTAU3

      

      REAL*4    RATIO(4),TMPPOS(3),XTMP,DUM
      INTEGER*4 NUMNEU(4)
      REAL*4    TMPRATIO(4)

      REAL*4    FNFIDMAS,RLU
      EXTERNAL  FNFIDMAS,RLU
      
      CHARACTER*8 CD,CT
      
C--Initialize FILE
      CALL KZINIT
      CALL SKOPENF(LUNO,1,'Z',IERR)
      IF (IERR.ne.0) THEN
         WRITE(*,*) "CAN not open FILE(LUN=",LUNO,")"
         STOP
      ENDIF

C--Initialize GEOMETRY
      CALL GEOSET

      DWALL=DWALLNE
#ifdef UPMU_RCK
*-----> upmu
      VOLKTON= FNFIDMAS(RROCK)  !--- rock volume
*     *** CALCULATE CYLINDER VOLUME
      CylVKTon = PI*((UPMU_DIAM**2)*RROCK-(RTKTK**2)*HITKTK/2)
     $     /1000./1000./1000.*dens_rck
      write(*,*)'init volume',volkton,cylvkton
      UPMU_EVFACT = CylVKton/VOLKTON    !-- event factor
      VOLKTON=VOLKTON*UPMU_EVFACT       !-- factored Volkton
*-----<
      WRITE(*,*) "VOLKTON FACT",VOLKTON,CylVKton,UPMU_EVFACT,RROCK
#else

      VOLKTON= FNFIDMAS(DWALL)
#endif
      IF (VOLKTON.LE.0.) THEN
         WRITE(*,*) "FIDUCAIAL MASS IS EQUAL OR LESS THAN 0..."
         STOP
      ENDIF

      ISOLR=ISOL
      SOLACTV=SOLACT
      CALL NEEVRATE(ISOLR,SOLACTV,EVRATE,RATIO)
C      NUMEVT=IFIX(EVRATE*YEAR*VOLKTON)

      NUMNEU(1)=IFIX(EVRATE*YEAR*VOLKTON*RATIO(1))
      NUMNEU(2)=IFIX(EVRATE*YEAR*VOLKTON*(RATIO(2)-RATIO(1)))
      NUMNEU(3)=IFIX(EVRATE*YEAR*VOLKTON*(RATIO(3)-RATIO(2)))
      NUMNEU(4)=IFIX(EVRATE*YEAR*VOLKTON*(1.-RATIO(3)))

      NUMEVT=NUMNEU(1)+NUMNEU(2)+NUMNEU(3)+NUMNEU(4)
      
      WRITE(6,*)   "---- PARAMETERS USED IN NEMAKE ---"
#ifdef UPMU_RCK
      WRITE(6,*) SOLACT,RROCK/100,VOLKTON,CYLVKTON,YEAR,NUMEVT,
     $             RATIO(1),NUMNEU(1),RATIO(2)-RATIO(1),
     $             NUMNEU(2),RATIO(3)-RATIO(2),NUMNEU(3),
     $             1.-RATIO(3),NUMNEU(4)

      WRITE(6,600) SOLACT,RROCK/100,VOLKTON,CYLVKTON,YEAR,NUMEVT,
     $             RATIO(1),NUMNEU(1),RATIO(2)-RATIO(1),
     $             NUMNEU(2),RATIO(3)-RATIO(2),NUMNEU(3),
     $             1.-RATIO(3),NUMNEU(4)

  600 FORMAT(' ',' SOL.ACTVITY     = ',F8.2,/,
     &          '  RROCK           = ',F10.2,/,
     &          '  KTON(FID.MAS.)  = ',F10.2,/,
     &          '  CYLVKTON(FID.MAS.)  = ',F10.2,/,
     &          '  YEAR    = ',F8.2 ,/,
     &          '  NUMBER OF EV    = ',I10 ,/,
     &          '  RATIO(neu_e)    = ',F8.3,' :# =',I8,/,
     &          '  RATIO(neu_ebar) = ',F8.3,' :# =',I8,/,
     &     '  RATIO(neu_mu)   = ',F8.3,' :# =',I8,/,
     &     '  RATIO(neu_mubar)= ',F8.3,' :# =',I8)

#else
      
      IF(ITAUFLG.EQ.1) THEN
         WRITE(6,*) "*** TAU MODE!!!! ***"
         WRITE(6,600) SOLACT,DWALL,VOLKTON,YEAR,NUMEVT,
     $        RATIO(1),NUMNEU(1),RATIO(2)-RATIO(1),
     $        NUMNEU(2),RATIO(3)-RATIO(2),NUMNEU(3),
     $        1.-RATIO(3),NUMNEU(4)
      ELSE
         WRITE(6,610) SOLACT,DWALL,VOLKTON,YEAR,NUMEVT,
     $        RATIO(1),NUMNEU(1),RATIO(2)-RATIO(1),
     $        NUMNEU(2),RATIO(3)-RATIO(2),NUMNEU(3),
     $        1.-RATIO(3),NUMNEU(4)
      ENDIF
  600 FORMAT(' ',' SOL.ACTVITY     = ',F8.2,/,
     &          '  WALL            = ',F8.3,/,
     &          '  KTON(FID.MAS.)  = ',F8.3,/,
     &          '  YEAR            = ',F8.3,/,
     &          '  NUMBER OF EV    = ',I8 ,/,
     &          '  RATIO(neu_e)    = ',F8.3,' :# =',I8,/,
     &          '  RATIO(neu_ebar) = ',F8.3,' :# =',I8,/,
     &     '  RATIO(neu_tau)   = ',F8.3,' :# =',I8,/,
     &     '  RATIO(neu_taubar)= ',F8.3,' :# =',I8)
  610 FORMAT(' ',' SOL.ACTVITY     = ',F8.2,/,
     &          '  WALL            = ',F8.3,/,
     &          '  KTON(FID.MAS.)  = ',F8.3,/,
     &          '  YEAR            = ',F8.3,/,
     &          '  NUMBER OF EV    = ',I8 ,/,
     &          '  RATIO(neu_e)    = ',F8.3,' :# =',I8,/,
     &          '  RATIO(neu_ebar) = ',F8.3,' :# =',I8,/,
     &     '  RATIO(neu_mu)   = ',F8.3,' :# =',I8,/,
     &     '  RATIO(neu_mubar)= ',F8.3,' :# =',I8)

#endif

#ifdef UPMU_RCK
      WRITE(6,*) "atrawetwet"
      NUMEVT=NUMNEU(3)+NUMNEU(4) !-- only numu event
      WRITE(*,*) 'TOTAL:',NUMEVT,'EVT'
#endif

      DO 100 I=1,NUMEVT
C         CALL RDMOUT(ISEED)
C         CALL RNRMBOUT(ISEED2)
         CALL RLUXAT(LUX,ISEED1,K1,K2)
         CALL RM48UT(ISEED2,K3,K4)
C         WRITE(6,*) ' SEED   : ',ISEED1,K1,K2,ISEED2,K3,K4
         CALL RMARUT(ISTAU1,ISTAU2,ISTAU3)
         CALL DATIMH(CD,CT)
#ifdef UPMU_RCK
         NUMREM     =NUMNEU(3)+NUMNEU(4) !-- only numu event
         TMPRATIO(3)=REAL(NUMNEU(3)                    )/REAL(NUMREM)
#else

         NUMREM     =NUMNEU(1)+NUMNEU(2)+NUMNEU(3)+NUMNEU(4)
         TMPRATIO(1)=REAL(NUMNEU(1)                    )/REAL(NUMREM)
         TMPRATIO(2)=REAL(NUMNEU(1)+NUMNEU(2)          )/REAL(NUMREM)
         TMPRATIO(3)=REAL(NUMNEU(1)+NUMNEU(2)+NUMNEU(3))/REAL(NUMREM)
#endif         
 10      CONTINUE
         IPAR=0
         
         XTMP=RLU(DUM)
#ifdef UPMU_RCK
         IF(                        XTMP.LE.TMPRATIO(3)) THEN
#else

         IF(                        XTMP.LE.TMPRATIO(1)) THEN
            IPAR= 12
            NUMNEU(1)=NUMNEU(1)-1
         ENDIF

         IF(XTMP.GT.TMPRATIO(1).AND.XTMP.LE.TMPRATIO(2)) THEN
            IPAR=-12
            NUMNEU(2)=NUMNEU(2)-1
         ENDIF

         IF(XTMP.GT.TMPRATIO(2).AND.XTMP.LE.TMPRATIO(3)) THEN
#endif 
           IPAR= 14
            IF(ITAUFLG.EQ.1) IPAR= 16
            NUMNEU(3)=NUMNEU(3)-1
         ENDIF

         IF(XTMP.GT.TMPRATIO(3)                        ) THEN
            IPAR=-14
            IF(ITAUFLG.EQ.1) IPAR=-16
            NUMNEU(4)=NUMNEU(4)-1
         ENDIF

         IF (IPAR.eq.0) THEN
            WRITE(6,*) "Warning!! IPAR was not defined(nemake.F)"
            GOTO 10
         ENDIF

#ifndef UPMU_RCK
         write(6,'(A8,I8,A2,A10,A2,A8)')
     $      'Event ',i,'  ',CD,'  ',CT
         write(6,'(A10,6I10)')
     $      '  Seeds:',iseed1,k1,k2,iseed2,k3,k4
#endif 
        IF (ITAUFLG.eq.1.and.abs(IPAR).ne.16) GOTO 100

#ifdef UPMU_RCK
         CALL NEMKEVNT(IPAR,UFLAG)
         IF(UFLAG.EQ.1) THEN
            STORE=STORE+1
            CALL KZWRIT(LUNO)
         ENDIF
c     Simulated volume of FC/PC is expanded to the entire water tank (2007,Dec. G.Mitsuka)
c#elif UPMU_WTR
c         CALL RNPOS_WTR(TMPPOS,DWALL)
c         CALL NEMKEVNT(IPAR,TMPPOS)
c         CALL KZWRIT(LUNO)
#else
         if (NEWATERVOL.eq.0) CALL RNPOS_WTRTNK(TMPPOS)    ! Full water volume
         if (NEWATERVOL.eq.1) CALL RNPOS(TMPPOS,DWALL)     ! ID-dwall
         if (NEWATERVOL.eq.2) CALL RNPOS_WTR(TMPPOS,DWALL) ! OD
         CALL NEMKEVNT(IPAR,TMPPOS)
         
         CALL KZWRIT(LUNO)
#endif
         CALL KZECLR

 100  CONTINUE
         
      CALL KZEND(LUNO)
      close(LUNO)

      RETURN
      END

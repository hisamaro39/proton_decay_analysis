      PROGRAM NEUT

c#include <f77/f77_floatingpoint.h>

#include "skhead.h"
#include "nework.h"
#include "nesolact.h"
#include "mcgenpar.h"      
#include "taucom.h"      
#include "necard.h"
#include "necardap.h"
#include "necardvc.h"
#include "neutmodel.h"


      REAL*4 POSI(3)
      REAL*4 DIR(3),GDIR(3)
      REAL*4 COST,PHI

#ifdef UPMU_RCK
       write(6,*) "********************"
       write(6,*) "FOR UPMU-ROCK"
       write(6,*) "********************"
#else
       write(6,*) "********************"
       write(6,*) "FOR FC,PC,UPMU-WATER"
       write(6,*) "********************"
#endif
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      LUNO = 20

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

c      CALL ieee_handler('set','common',SIGFPE_ABORT)
c      CALL ieee_handler('set','division',SIGFPE_ABORT)
c      CALL ieee_handler('set','invalid',SIGFPE_ABORT)

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

C -- this will be defined in the necardvc.F
C      SK_GEOMETRY = 2
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC      

      ! load configuration for interactions
      CALL NECARD
      ! load configuration for atmospheric neutrino
      CALL NECARDAP
      CALL NECARDVC

      ITAUFLGCORE = ITAUFLG

C--Fill interaction models
      call nefillmodel

C--Fill liblary version
      call nefillver

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC            

      LUNO=20
      YEAR = NEYEARS
      CALL NEMAKE(LUNO,YEAR)

 10   CONTINUE

C      CALL RDMOUT(ISEED)
C      CALL RNRMBOUT(ISEED2)
      CALL RLUXAT(LUX,ISEED1,K1,K2)
      WRITE(6,*) ' RANDOM NUMBER SEED '
      WRITE(6,*) ' RANLUX : ',LUX,ISEED1,K1,K2
      CALL RM48UT(ISEED2,K1,K2)
      WRITE(6,*) ' RM48 : ',ISEED2,K1,K2
      CALL RMARUT(ISTAU1,ISTAU2,ISTAU3)

      CALL NECHADIR(COST,PHI,DIR,GDIR)

      STOP
      END

***********************************************************************     
*     --------------------------
      SUBROUTINE RNPOS_WTRTNK(POS)
*     --------------------------
*     
*     ( purpose )
*       SET RANDOM POSITION
*
*     ( input ) 
*       None
*
*     ( output ) 
*       POS(3) ( cm )
*
*     ( creation date and author )
*       2007.12.05 ; First creation based on rnpos.F by G.Mitsuka
*
***********************************************************************     
      implicit none

#include "geotnk.h"

      Real*4 POS(3)
c      DIMENSION POS(3)
      Real*4 RLU, SSS
      external RLU
      INTEGER IDUM

*     Set z position(From the bottom to the top of water tank)
 100  POS(3)=HITKTK*(2.*RLU(IDUM)-1.)
      if(abs(POS(3)).gt.ZPTKTK) goto 100

*     Set x and y position(To the wall of water tank)
 101  POS(1)=(2.*RLU(IDUM)-1.)*RTKTK
      POS(2)=(2.*RLU(IDUM)-1.)*RTKTK
      SSS=SQRT(POS(1)**2+POS(2)**2)
      IF(SSS.GT.RTKTK)GO TO 101

      RETURN
      END

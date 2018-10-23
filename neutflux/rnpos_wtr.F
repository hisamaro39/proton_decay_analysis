***********************************************************************     
*     --------------------------
      SUBROUTINE RNPOS_WTR(POS,dwall)
*     --------------------------
*     
*     ( purpose )
*       SET RANDOM POSITION
*
*     ( output ) 
*       POS(3) ( cm )
*
*     ( creation date and author )
*
***********************************************************************     
#include "geotnk.h"
      DIMENSION POS(3)

100   POS(1)=(2.*RLU(DUM)-1.)*RTKTK
      POS(2)=(2.*RLU(DUM)-1.)*RTKTK
      POS(3)=(2.*RLU(DUM)-1.)*ZPTKTK

      SS=SQRT(POS(1)**2+POS(2)**2)
c      IF(abs(pos(3)).lt.ZPINTK) GO TO 100
c      if(ss.gt.RTKTK .or. ss.lt.RINTK) goto 100
      if(ss.gt.RTKTK) goto 100 
      IF(abs(pos(3)).lt.ZPINTK .and. ss.lt.RINTK) GO TO 100

      RETURN
      END

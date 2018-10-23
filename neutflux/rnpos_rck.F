***********************************************************************     
*     --------------------------
      SUBROUTINE RNPOS_RCK(POS,radrck,theta,phi)
*     --------------------------
*     
*     ( purpose )
*       SET RANDOM POSITION
*
*     ( output ) 
*       POS(3) ( cm )
*
*     ( creation date and author )
*     radius(cm)
***********************************************************************     
      implicit none
#include "geotnk.h"
#include "upmu.h"
      real pos(3),sss,ss
      real radrck,dum
      real RLU
      EXTERNAL RLU

      real ipos(3),trans(3)
      real theta,phi,r

*---   IPOS(3): initial position
 200  ipos(1)=(2.*rlu(dum)-1)*UPMU_DIAM
      ipos(2)=(2.*rlu(dum)-1)*UPMU_DIAM
      ipos(3)=rlu(dum)*radrck
      r=sqrt(ipos(1)**2+ipos(2)**2)
      if(r.gt.UPMU_DIAM) goto 200
**---  rotate y axis(theta)
      trans(1)=ipos(1)*cos(theta)-ipos(3)*sin(theta)
      trans(2)=ipos(2)
      trans(3)=ipos(1)*sin(theta)+ipos(3)*cos(theta)
*---  rotate z axis(phi)
      pos(1)=trans(1)*cos(phi)-trans(2)*sin(phi)
      pos(2)=trans(1)*sin(phi)+trans(2)*cos(phi)
      pos(3)=trans(3)
*---  inner tank & out of range check(reject)
      SSS=SQRT(POS(1)**2+POS(2)**2+pos(3)**2)
      IF(SSS.GT.radrck) GO TO 200

      SS=SQRT(POS(1)**2+POS(2)**2)
      IF(SS.lt.rtktk .and. abs(pos(3)).lt.ZPtkTK) GO TO 200

      return
      end


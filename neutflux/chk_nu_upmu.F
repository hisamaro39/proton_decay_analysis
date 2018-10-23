*********************************************************
      Subroutine chk_nu_upmu( pos, energy, dir, IWRT_UPMU )
*
*     (PURPOSE)  CHECK NU THROUGH NEAR THE DETECTOR
*     (INPUT)
*        POS:    NU  INTERACTION POSITION
*        ENERGY: NU  ENERGY(MeV)
*        DIR:    NU  DIRECTION
*     (OUTPOT)
*        IWRT_UPMU:  NEAR(=1)OR FAR(=0)
*     (EDITOR)
*        C.SAJI
********************************************************  
      Implicit none
#include "upmu.h"
      Real pos(3),dir(3),pdir(3)
      real energy,distance,theta,theta_limit,range,tmp_acos
      real ang_mu,ang_nu,X0
      real*8 range1,E0
      Real Scal
      EXTERNAL SCAL
      integer IWRT_UPMU

      X0=10.5                     !-- mu multiple scatter param.
      E0=1.6                    !-- MU ENERGY THRESHOLD
      CALL V3MULT(-1.,POS,PDIR)   !-- DIR. OF DETECTOR
      distance =sqrt( Scal( pos,pos ))    !-- DISTANCE
      CALL hrange(dble(energy/1000.),E0,range1) !-- RANGE
      range=range1/dens_rck*3

C      write(11,*) 'chk',energy,range

      If( distance.le.UPMU_RMIN ) then
         iwrt_upmu = 1              !-- int. pos. is inner of UPMU_RMIN
      Elseif( distance.le.(range+UPMU_RMIN) ) then
         tmp_acos=Scal(dir,pdir)/SQRT(SCAL(DIR,DIR))/distance 
         if( tmp_acos.GE.1. ) then
            theta = 0.0
         else
            theta = Acos(tmp_acos)
         endif
C         do i=1,3
C            pd(i)=pdir(i)/sqrt(Scal(pdir,pdir))
C            di(i)=dir(i)/SQRT(SCAL(DIR,DIR))
C         enddo
         ang_mu = 13.6/energy*sqrt(distance/x0)*
     &        (1.+0.038*log(distance/x0))*sqrt(2.) !-- mu scatter
         ang_nu = 2.6*sqrt(100000./energy)/180.*PI*sqrt(2.) !-- nu

C         write(1,*) energy,ang_nu*180/PI,ang_mu*180/PI,
C     &        1*(ang_mu+ang_nu)*180/PI
C         write(2,*) 'DM',UPMU_RMIN,distance,ang_nu,ang_mu

***********************************************************************
*---       *** mu dir limit(Multiple scattering through smaill angles)
*---       *** 6*sigma
           theta_limit = Asin( UPMU_RMIN/distance )+6*(ang_mu+ang_nu)
***********************************************************************
C           write(2,*) 'E',energy,theta, theta_limit, 
C     &      Asin( UPMU_RMIN/distance )

           If( theta.le.theta_limit ) then
              iwrt_upmu = 1       !-- NU DIR. OK
           else
              iwrt_upmu = 0       !-- NU DIR. BAD SKIP...
           Endif
        Else
           iwrt_upmu = 0          !-- NU SHORT RANGE SKIP...
        Endif
      Return
      End

      Real Function Scal( vec1, vec2 )
      Implicit none
      Real vec1(3),vec2(3)
      Scal = vec1(1)*vec2(1)+vec1(2)*vec2(2)+vec1(3)*vec2(3)
      Return
      End


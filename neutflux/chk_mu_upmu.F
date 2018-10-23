*********************************************************
      Subroutine chk_mu_upmu(pos,pm,dir,iwrt1_upmu)
*
*     (PURPOSE)  CHECK MU THROUGH NEAR THE DETECTOR
*     (INPUT)
*        POS:    MU  INTERACTION POSITION(cm)
*        ENERGY: MU  MOMENTUM(MeV)
*        DIR:    MU  DIRECTION
*     (OUTPOT)
*        IWRT_UPMU1:  NEAR(=1)OR FAR(=0)
*     (EDITOR)
*        C.SAJI
********************************************************  
      Implicit none
#include "upmu.h"
      INTEGER IWRT1_UPMU
      Real pos(3),energy,dir(3),pdir(3)
      Real distance,range,theta,theta_limit,ang_sigma,tmp_acos
      Real X0,pm
      Real Scal
      EXTERNAL SCAL
      real*8 range1,E0

      X0=10.5                   !-- FOR MU MULTIPLE SCATTER
      E0=1.6                  !-- MU ENERGY THRESHOLD
      energy =sqrt(pm**2+105**2)
      CALL V3MULT(-1.,POS,PDIR) !--DIR. OF DETECTOR
      distance =sqrt(Scal(pos,pos )) !-- DISTANCE
      CALL hrange(dble(energy/1000.),E0,range1) !-- RANGE
      range=range1/dens_rck*3.

      If( distance.le.UPMU_RMIN ) then
         iwrt1_upmu = 1         !-- INT. POS IS NEAR BY TANK
      Elseif( distance.le.(range+UPMU_RMIN) ) then
         tmp_acos = Scal(dir,pdir)/sqrt(Scal(dir,dir))/distance 
         if( tmp_acos.GE.1. ) then
            theta = 0.0
         else
            theta = Acos(tmp_acos)
         endif

C         do i=1,3
C            pd(i)=pdir(i)/sqrt(Scal(pdir,pdir))
C            di(i)=dir(i)/sqrt(Scal(dir,dir))
C         enddo
*---       *** mu dir limit(Multiple scattering through smaill angles)
*---       *** 4*sigma
         ang_sigma=13.6/PM*sqrt(distance/x0)* 
     &        (1.+0.038*log(distance/x0))*sqrt(2.)
         theta_limit = Asin( UPMU_RMIN/distance )+8*ang_sigma
         If( theta.le.theta_limit ) then
              iwrt1_upmu = 1
           else
              iwrt1_upmu = 0
C     write(*,*) ' neutrino direction bad. Skip... '
           Endif
        Else
           iwrt1_upmu = 0
C     write(*,*) ' neutrino short range. Skip... '
        Endif
      Return
      End




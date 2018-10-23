***********************************************************************     
*     --------------------------
      SUBROUTINE dirtrans(dir,theta,phi)
*     --------------------------
*     
*     ( purpose )
*       transter the direction from x,y,z to cosz,theta
*
*     ( output ) 
*       POS(3) ( cm )
*
*     ( creation date and author )
*
***********************************************************************     
      implicit none

      real*4 dir(3),dirt(3),dirtt(3),theta,phi
      real pi
      real*4 costh,cosphi
      integer i
      
      DATA PI/3.141593/      

      do i=1,3
         dirt(i)=dir(i)
      enddo

*     **** THETA
      costh=-dirt(3)
      theta=acos(costh)
*     **** PHI
      if((dirt(1)**2+dirt(2)**2).LE.0.) then
         phi=0.
      else
         cosphi=abs(dirt(1))/sqrt(dirt(1)**2+dirt(2)**2)
         phi=acos(cosphi)
      endif

      if(dirt(1).gt.0.and.dirt(2).gt.0) then 
         phi=phi
      elseif(dirt(1).gt.0.and.dirt(2).le.0) then 
         phi=-phi
      elseif(dirt(1).le.0.and.dirt(2).gt.0) then 
         phi=PI-phi
      elseif(dirt(1).le.0.and.dirt(2).le.0) then 
         phi=-(PI-phi)
      endif

      dirtt(1)=cos(phi)*sin(theta)
      dirtt(2)=sin(phi)*sin(theta)
      dirtt(3)=-cos(theta)

      RETURN
      END





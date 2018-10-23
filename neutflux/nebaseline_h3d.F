      subroutine nebaseline_h3d(path,nutype,costheta,E_neu,d)
      
      implicit none

CC output
      real path(20)

C production height
      real hgt(20)
      real P
      real costheta,E_neu,d
      integer nutype

* Return random baseline at given costheta and neutrino energy,
*  for a detector which is a distance d in km above the surface
* of the earth.
*
* Calls "neprodhgt_h3d.F" to get the slant height in the atmosphere.
*
* 1-electron neutrino
* 2-muon neutrino
* 3-anti electron neutrino
* 4-anti muon neutrino

      real Q ! Distance to surface along zenith angle 
      real cosalpha ! Angle at surface
      
      real a,b,c,x
      integer ii
      
      real R
      parameter (R=6371.000) ! Radius of the Earth in km 

* Define values for quadratic solution 
      a = 1
      x = R+d
  
      if (costheta.lt.0.) then
         b = 2.*x*costheta
      else
         b = -2.*x*costheta
      end if
      c = x*(x)-R*R
  
      if (b*b.lt.4.*a*c) then
* Only true very near horizon.  Use tangent value there.
         Q = sqrt(2.*d*x)
      else
         if (costheta.gt.0) then
            Q = (-b - sqrt(b*b-4.*a*c))/(2.*a)
         else
            Q = (-b + sqrt(b*b-4.*a*c))/(2.*a)
         end if
      end if
      cosalpha = abs((Q*Q+R*R-x*x)/(2.*Q*R))
      call neprodhgt_h3d(hgt,nutype,costheta,E_neu)
 
* Now return the baseline

      if(abs(costheta).le.sqrt((x**2-R**2)/x**2)) then
         do ii=1,20
            path(ii) = -x*costheta + sqrt((x*costheta)**2+(R+hgt(ii))**2-x**2)
         end do
      else
         do ii=1,20
            P = sqrt((R*cosalpha)**2+(R+hgt(ii))**2-R**2)-R*cosalpha
            if (costheta.lt.0.) then
               path(ii) = Q+P
            else 
               path(ii) = P-Q
            end if
         end do
      end if

      return 
      end




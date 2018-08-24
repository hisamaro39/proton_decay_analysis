      subroutine neprodhgt_h3d(hgt,nutype, costheta, E_nu)

      implicit none

      real h3d(19,31,20,4)
      save h3d
C  h3d: probabiliy, E_nu, zenith angle (cos(theta)), particle ID
c  probability;  5%   10%  11% ..... 90%  95%
c  zenith angle; 01 = -1.0 < cosz < -0.9, 02 = -0.9 < cosz < -0.8,  ....
c                19 =  0.8 < cosz <  0.9, 20 =  0.9 < cosz <  1.0, 
C  particle ID; (numu,numu_bar,nue.nue_bar)=(1,2,3,4) (or =(1,2,-1,-2))

CC TABLE
      character*60 table
      data table/'/home/atmpd/const/honda3d_prodhgt'/
      character*256 filename
      logical lexist
      integer jj, iflavor, izenith, ienergy
      real edum
      integer icallprod
      data icallprod/0/
CC Output
      real hgt(20)
CC Input
      integer nutype
      real costheta, E_nu
CC Local
      integer ii
      integer costag,enetag,nutag
      real delta_h
      real dhe,delta_e,de
      real dhc,delta_cos
      
C --  read table
      if(icallprod .eq. 0) then
         icallprod = 1
         inquire(file=table,exist=lexist)
         if ( lexist ) then
            open(34,file=table,status='old')
         else
            call findconsts(filename,'production height',
     $           'honda3d_prodhgt',0)
            open(34,file=filename,status='old')
         endif
         do iflavor=1,4
            do izenith=1,20
               do ienergy=1,31                  
                  read(34,*) edum
     &                 ,(h3d(jj,ienergy,izenith,iflavor),jj=1,19)
               end do
            end do
         end do
         close(34)
      end if

      costag = 1
      costag = int((costheta+1.0)*10.)+1
      if(costag.gt.20) costag=20 
      delta_cos = costheta - ((costag-11.0)/10.0+0.05)

      enetag = int(10.0*log10(E_nu) + 11.0)
      if(enetag .le.  0) enetag =  1
      if(enetag .gt. 31) enetag = 31
      de = 10**(((enetag+1)-11.0)/10.0) - 10**((enetag-11.0)/10.0)
      delta_e =  E_nu - 10**((enetag-11.0)/10.0)

      if(nutype.eq.1) nutag = 3
      if(nutype.eq.2) nutag = 1
      if(nutype.eq.3 .or. nutype.eq.-1) nutag = 4
      if(nutype.eq.4 .or. nutype.eq.-2) nutag = 2

      do ii=1,19

         delta_h = 0.0
         
         if(enetag.lt.31 .and. delta_e.gt.0.0) then
            dhe = h3d(ii,enetag+1,costag,nutag)
     &           -h3d(ii,enetag,costag,nutag)
            delta_h = delta_h + delta_e * dhe/de
         end if
         
         if((costag.eq.1.and.delta_cos.le.0.)
     &        .or.(costag.eq.20.and.delta_cos.ge.0.)) then
            dhc = 0.0
         else
            if(delta_cos.ge.0.) then
               dhc = h3d(ii,enetag,costag+1,nutag)
     &              -h3d(ii,enetag,costag,nutag)
            else
               dhc = h3d(ii,enetag,costag-1,nutag)
     &              -h3d(ii,enetag,costag,nutag)
            end if
            delta_h = delta_h + abs(delta_cos) * dhc/0.1
         end if

         hgt(ii) = h3d(ii,enetag,costag,nutag) + delta_h

      end do

      hgt(20) = hgt(19)

      return
      end


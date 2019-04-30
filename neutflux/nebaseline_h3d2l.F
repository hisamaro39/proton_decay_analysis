      subroutine nebaseline_h3d2l(path,nutype, costheta, E_nu)

      implicit none

      real h3d2l(16,21,40,4)
      save h3d2l
C  h3d2l: probabiliy, E_nu, zenith angle (cos(theta)), particle ID
c  energy of neutrino; 0.1GeV ... 100GeV, 5bins in a decade, 16bins in total 
c    comment; 12th,13th,14th,15th and 16th bin of nue has path length of 0km
c                         -> not used
c  probability;  1% 5%   10%  15% ..... 90%  95% 99%, in total 21bins
c  zenith angle; 01 =  1.0  < cosz < 0.95
c                02 =  0.95 < cosz < 0.90
c                 ...
c                40 = -0.95 < cosz < -1.00 , in total 40bins
C  particle ID; (numu,numu_bar,nue,nue_bar)=(1,2,3,4)
c                                    input =(14,-14,12,-12)

CC TABLE
      character*60 table
      data table/'neutflux/dist-low-0.90-aa.d'/
      integer jj, iflavor, izenith, iprob
      character*500 chdum
      real pdum
      integer icallprod
      data icallprod/0/
CC Output (km)
      real path(20)
CC Input (PDG code, cos(thenith), GeV)
      integer nutype
      real costheta, E_nu
CC Local
      integer ii
      integer costag,enetag,nutag
      real delta_h
      real dhe,delta_e
      real dhc,delta_cos
      
C --  read table
      if(icallprod .eq. 0) then
         icallprod = 1
         open(34,file=table,status='old')
         do iflavor=1,4
            do izenith=1,40
               read(34,*) chdum
               read(34,*) chdum
               do iprob=1,21                  
                  read(34,*) pdum
     &                 ,(h3d2l(jj,iprob,41-izenith,iflavor),jj=1,16)
               end do
            end do
         end do
         close(34)
      end if

      costag = 1
      costag = int((costheta+1.0-0.025)*20.)+1
      if(costag.gt.40) costag=40 
      delta_cos = (costheta+1.0-0.025)*20.+1.0 - float(costag)

      enetag = int(5.0*log10(E_nu)+6)
      if(enetag .le.  0) enetag =  1
      if(enetag .gt. 11) enetag = 11
      delta_e =  5.0*log10(E_nu)+6 - float(enetag)

      if(nutype.eq.14) nutag = 1
      if(nutype.eq.-14) nutag = 2
      if(nutype.eq.12) nutag = 3
      if(nutype.eq.-12) nutag = 4

      do ii=1,20

         delta_h = 0.0
         
         if(enetag.lt.11 .and. delta_e.gt.0.0) then
            dhe = (h3d2l(enetag+1,ii+1,costag,nutag)
     &           +h3d2l(enetag+1,ii,costag,nutag))/2.0
     &           -(h3d2l(enetag,ii+1,costag,nutag)
     &           +h3d2l(enetag,ii,costag,nutag))/2.0
            delta_h = delta_h + delta_e * dhe
         end if
         
         if(costag.lt.40 .and. delta_cos.ge.0.) then
            dhc = (h3d2l(enetag,ii+1,costag+1,nutag)
     &           +h3d2l(enetag,ii,costag+1,nutag))/2.0
     &           -(h3d2l(enetag,ii+1,costag,nutag)
     &           +h3d2l(enetag,ii,costag,nutag))/2.0
            delta_h = delta_h + delta_cos * dhc
         end if

         path(ii) = (h3d2l(enetag,ii+1,costag,nutag)
     &        +h3d2l(enetag,ii,costag,nutag))/2. + delta_h
         path(ii) = path(ii)/1000.

      end do

      return
      end


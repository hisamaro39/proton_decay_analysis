
      SUBROUTINE NERDNETARG

      IMPLICIT NONE

#include "skhead.h"
#include "nework.h"
#include "necard.h"
#include "nrcard.h"
#include "nefillver.h"
#include "neutmodel.h"
#include "neutparams.h"
#include "neutcrs.h"
#include "vcwork.h"
#include "posinnuc.h"
#include "rscons.h"

      integer*4 NSEGM
      PARAMETER (NSEGM=8)

      integer*4 LSEGMTBL(10)
      DATA LSEGMTBL/ 3, 6,11, 7, 5,27,27, 10, 0, 0/

      integer*4 IDATA(500)
      real*4    RDATA(500)

      EQUIVALENCE(IDATA(1),RDATA(1))

      CHARACTER*32 CNAME

      integer*4 IERR,ISEGM,NDATA,I

      CNAME='NEUTINFO'

      ierr = 0
      CALL KZBLOC(CNAME,IERR)
      if (ierr.eq.0) then 
         goto 999
      endif

C--------------------------      
C- EMB
C-    + 1 : Length of this segment
C     + 2 : Number of segments
C-    + 3 : Length of 1st segment  : Versions of libraries
C     + 4 : Length of 2nd segment  : Flux & Geometry
C     + 5 : Length of 3rd segment  : Interaction and nuclear effects
C     + 6 : Length of 4th segment  : Reaction models
C     + 7 : Length of 5th segment  : Target information
C     + 8 : Length of 6th segment  : Cross-section weights ( neutrino )
C     + 9 : Length of 7th segment  : Cross-section weights ( anti neutrino )
C     +10 : Length of 8th segmetn  : Total cross-section
C

      ISEGM   = 0
      call kzget0(CNAME,ISEGM,NDATA,IDATA)

      DO 10 I=1,NSEGM
         if (IDATA(2+I).ne.LSEGMTBL(I)) then
            write (*,*) CNAME,':Bank length for seg #',i,
     $           ' (',IDATA(2+I),') is different from expected:',
     $           ' (',LSEGMTBL(I),').'
            stop
         endif
 10   continue

C--------------------------      
C- Segment 1
C    Version of libraries
C     + 1 : NeutCore
C     + 2 : Nuceff
C     + 3 : Nucleon rescattering
C     

      ISEGM = 1
      
      call kzget0(CNAME,ISEGM,NDATA,IDATA)
      if (NDATA.ne.LSEGMTBL(ISEGM)) then
         write (*,*) CNAME,':Bank length for seg #',ISEGM,
     $        ' (',NDATA,') is different from expected',
     $        ' (',LSEGMTBL(ISEGM),').'
         stop
      endif

      corev = IDATA(1)
      nucev = IDATA(2)
      nuccv = IDATA(3) 

C--------------------------      
C- Segment 2
C    Flux & Geometry
C     + 1 : SK Geometry
C     + 2 : Solar Activity on/off ( for ATMPD )
C     + 3 : Solar Activity param. ( for ATMPD )
C     + 4 : Neut Flux type        ( for ATMPD )
C     + 5 : Years                 ( for ATMPD )
C     + 6 : ITAUFLG               ( for ATMPD )

      ISEGM    = 2

      call kzget0(CNAME,ISEGM,NDATA,IDATA)
      if (NDATA.ne.LSEGMTBL(ISEGM)) then
         goto 9999
      endif

      SK_GEOMETRY = IDATA(1)

      
C--------------------------      
C- Segment 3
C    Interaction
C     + 1 : Fermi momentum flag
C     + 2 : Pauli blocking flag
C     + 3 : Nuclear effect flag
C     + 4 : Interaction mode
C     + 5 : Fermi surface momentum
C     + 6 : Maximum Fermi momentum
C     + 7 : Initial nuclear potential        
C     + 8 : Final nuclear potential
C     + 9 : Formation zone flag
C     +10 : Nucleon rescattering flag     
C     +11 : Factor for nucleon rescattering prob.
C      

      ISEGM    = 3

      call kzget0(CNAME,ISEGM,NDATA,IDATA)
      if (NDATA.ne.LSEGMTBL(ISEGM)) then
         goto 9999
      endif

      NEFRMFLG = IDATA(1) 
      NEPAUFLG = IDATA(2) 
      NENEFO16 = IDATA(3) 
      NEMODFLG = IDATA(4) 
      PFSURF   = RDATA(5) 
      PFMAX    = RDATA(6) 
      VNUINI   = RDATA(7) 
      VNUFIN   = RDATA(8) 
      IFORMLEN = IDATA(9) 
      NUCRESCAT= IDATA(10)
      XNUCFACT = RDATA(11)

C--------------------------     
C- Segment 4
C       Reaction models
C     + 1 : QE model
C     + 2 : Single pi model
C     + 3 : Coherent pi model
C     + 4 : Deep Inelastic scattering model
C     + 5 : QE potential model
C     + 6 : QE Axial mass
C     + 7 : Single pi Axial mass
C

      ISEGM    = 4

      call kzget0(CNAME,ISEGM,NDATA,IDATA)
      if (NDATA.ne.LSEGMTBL(ISEGM)) then
         goto 9999
      endif

      MDLQE   = IDATA(1)
      MDLSPI  = IDATA(2)
      MDLCOH  = IDATA(3)
      MDLDIS  = IDATA(4)
      MDLQEAF = IDATA(5)     
      XMAQE   = RDATA(6)
      XMASPI  = RDATA(7)

      MODELDIS = MDLDIS
      MODELCOH = MDLCOH

C--------------------------      
C- Segment 5
C       Target information
C     + 1 : # of bound neutron
C     + 2 : # of bound proton
C     + 3 : # of free proton
C     + 4 : # of nucleons in an atom

      ISEGM    = 5

      call kzget0(CNAME,ISEGM,NDATA,IDATA)
      if (NDATA.ne.LSEGMTBL(ISEGM)) then
         goto 9999
      endif

      NUMBNDN = IDATA(1)
      NUMBNDP = IDATA(2)
      NUMFREP = IDATA(3)
      NUMATOM = IDATA(4)
      IBOUND  = IDATA(5)

C--------------------------      
C- Segment 6
C       Cross-section weights ( neutrino )
C     + 1 : 
C      |  : Cross-section weight for neutrino
C     +27 : 

      ISEGM    = 6

      call kzget0(CNAME,ISEGM,NDATA,IDATA)
      if (NDATA.ne.LSEGMTBL(ISEGM)) then
         goto 9999
      endif

      DO 30 I=1,27
         CRSNEUT(I) = RDATA(I)
 30   continue

C--------------------------      
C- Segment 7
C       Cross-section weights ( anti neutrino )
C     + 1 : 
C      |  : Cross-section weight for anti neutrino
C     +27 : 
C

      ISEGM    = 7

      call kzget0(CNAME,ISEGM,NDATA,IDATA)
      if (NDATA.ne.LSEGMTBL(ISEGM)) then
         goto 9999
      endif

      DO 40 I=1,27
         CRSNEUTB(I) = RDATA(I)
 40   continue

C--------------------------      
C- Segment 8
C       Total cross-section
C     + 1 : Energy of neutrino  ( GeV )
C     + 2 : Cross-section       (^10-38cm^2)


      ISEGM    = 8

      call kzget0(CNAME,ISEGM,NDATA,IDATA)
      if (NDATA.ne.LSEGMTBL(ISEGM)) then
         goto 9999
      endif

      CRSENERGY = RDATA(1)
      TOTCRSNE  = RDATA(2)
      
      DO 50 I=3,NDATA
         DIFCRSNE(I-2) = RDATA(I)
 50   continue

      return
      
 999  continue
      if (ierr.ne.0) then
         write(*,*) 'Failed to create BANK : ',CNAME
      endif

      return

 9999 continue
      write (*,*) CNAME,':Bank length for seg #',ISEGM,
     $     ' (',NDATA,') is different from expected',
     $     ' (',LSEGMTBL(ISEGM),').'
      stop
      
      end

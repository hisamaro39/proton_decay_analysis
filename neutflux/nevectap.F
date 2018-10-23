**********************************************************************
*     ------------------------------
      SUBROUTINE NEVECTAP(IP,POS,PMOM,IERR)
*     ------------------------------
*     
*     ( purpose )
*       VECTOR GENERATION ( NEUTRINO INTERACTION )
*
*     ( input ) 
*       IP      : PARICLE CODE
*       POS(3)  : VERTEX POSITION
*       PMOM(3) : MOMENTUM ( MEV/C)
*
*     ( output ) 
*       COMMON VCWORK,VCVRTX
*
*     ( creation date and author )
*       1995.02.25 ; First version by K.Kaneyuki
*       1996.01.11 ; K.Kaneyuki  NENEFO16 = 1 -> no nuclear effect in O16
*       1998.03.07 ; add nek0(k0 -> KL or KS) by J.K.
*       2004.07.06 ; save interaction position in O16 by M.Ishitsuka
*       2007.08.22 ; G.Mitsuka adds NEMKMODELBK
*       2009.12.09 ; R.Tacik, P. de Perio
*                    - Add nucleon ejection by pi/Delta absorption
*
**********************************************************************
      implicit none

#include "mcgenpar.h"
#include "necardap.h"
#include "necard.h"
#include "fsihist.h"

      integer*4 IP
      integer*4 IERR

      REAL POS(3),PMOM(3)
      INTEGER ICALL/0/

      INTEGER*4 LUX,ISEED1,K1,K2
      INTEGER*4 ISEED2,K3,K4

      real*4    EVPIPROB
      external  EVPIPROB
C
C
      ierr = 0
      write(*,*) '-----------------------------------------------------'

C
C --  CLEAR COMMON
C
      CALL VCCLCM
      CALL NECLPOS
C     
C --  MAKE HEADER
C
      CALL MCMKHD
      CALL MCMKMH
      CALL MCMKWT
C     
C --  MAKE VECTOR
C      
      CALL NEVENT(IP,PMOM,IERR)
      if (ierr.ne.0) then
         call KZBDEL('HEADER')
         call KZBDEL('MCHEAD')
         call KZBDEL('MCPARMWATR')
         return
      endif
C
C --  MAKE NEUT BANK
C
      CALL NEMKNEBK(POS)
      CALL NEMKMODELBK
      CALL NEMKCRSBK
C
C --  COPY PARAMETER FROM NEWOK TO VCWORK, VCVRTX
C
      CALL NECPNEWK(POS)
C
C --  TRACE IN O16
C
C      IF (NUCFLG1.EQ.0) CALL NEFFECT

      CALL RLUXAT(LUX,ISEED1,K1,K2)
      CALL RM48UT(ISEED2,K3,K4)
      WRITE(6,*) ' SEED (Before NEFFECT): ',ISEED1,K1,K2,ISEED2,K3,K4

      CALL EFCLFSI
      IF (NENEFO16.NE.0) then
         write(*,*) 'STOP!'
      endif
      IF (NENEFO16.EQ.0) CALL NEFFECT

      CALL RLUXAT(LUX,ISEED1,K1,K2)
      CALL RM48UT(ISEED2,K3,K4)
      WRITE(6,*) ' SEED (After NEFFECT): ',ISEED1,K1,K2,ISEED2,K3,K4
C
C
C --  CONSIDER NULCEAR EXCITATION OF GAMMA
C
C      CALL NEMKGAM
      if(NUMATOM.eq.16) CALL NEMKGAM ! only 16O

      CALL RLUXAT(LUX,ISEED1,K1,K2)
      CALL RM48UT(ISEED2,K3,K4)
      WRITE(6,*) ' SEED (After NEMKGAM): ',ISEED1,K1,K2,ISEED2,K3,K4

C
C --  decay omega meson 
C
      CALL NEOMEGA

      CALL RLUXAT(LUX,ISEED1,K1,K2)
      CALL RM48UT(ISEED2,K3,K4)
      WRITE(6,*) ' SEED (After NEOMEGA): ',ISEED1,K1,K2,ISEED2,K3,K4

C
C --  nucleon rescattering
C
      CALL RLUXAT(LUX,ISEED1,K1,K2)
      CALL RM48UT(ISEED2,K3,K4)
      WRITE(6,*) ' SEED (Before NRINTR): ',ISEED1,K1,K2,ISEED2,K3,K4

      call nrintr

      CALL RLUXAT(LUX,ISEED1,K1,K2)
      CALL RM48UT(ISEED2,K3,K4)
      WRITE(6,*) ' SEED (After NRINTR): ',ISEED1,K1,K2,ISEED2,K3,K4
C      
C --  nucleon ejection following pion/delta absorption
C
      call absneject
C
C --  radiative correction
C     
      call radcorr

C --  EXECUTE TAU DECAY (needs to be after NEMKNPOS)
C
      CALL TAUVECT

C
C --  translate K0 -> K LONG or K SHORT 
C
      CALL NEK0

      CALL RLUXAT(LUX,ISEED1,K1,K2)
      CALL RM48UT(ISEED2,K3,K4)
      WRITE(6,*) ' SEED (After NEK0): ',ISEED1,K1,K2,ISEED2,K3,K4

C
C -- Calculate FSIPROB after all possible instances of EFTRACE
C
      FSIPROB=EVPIPROB()
C     
C     
C --  MAKE BANK FOR OUTPUT
C
C
C --- MAKE NUCPOS BANK
      CALL NEMKNPOS

      CALL VCMKVC
      CALL VCMKVX
      CALL NEMKFSIBK
      CALL NEMKNETARG

C      
      RETURN
      END

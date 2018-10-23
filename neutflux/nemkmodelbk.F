***********************************************************************
*     ------------------------
      SUBROUTINE NEMKMODELBK
*     ------------------------
*     
*     ( purpose )
*       Create NEUT bank for reaction models etc...
*
*     ( input ) 
*       COMMON block
*
*     ( output )
*       ZBS bank
*
*     ( creation date and author )
*       2007.08.22 G.Mitsuka 
*                  makes new bank for flag of reaction models
*
***********************************************************************

#include "necard.h"
#include "skhead.h"
#include "nrcard.h"
#include "neutparams.h"
#include "neutmodel.h"
#include "nefillver.h"

      ! neut version
      INTEGER IDATAN(3)
      REAL    RDATAN(3)
      CHARACTER*32 CNAMEN
      EQUIVALENCE (IDATAN(1),RDATAN(1))

      ! flux and geometry
      INTEGER IDATAF(6)
      REAL    RDATAF(6)
      CHARACTER*32 CNAMEF
      EQUIVALENCE (IDATAF(1),RDATAF(1))

      ! nuclear effects
      INTEGER IDATAC(11)
      REAL    RDATAC(11)
      CHARACTER*32 CNAMEC
      EQUIVALENCE (IDATAC(1),RDATAC(1))

      ! reaction models
      INTEGER IDATAM(10)
      REAL    RDATAM(10)
      CHARACTER*32 CNAMEM
      EQUIVALENCE (IDATAM(1),RDATAM(1))


C--- NEUT version
      IDATAN(1)=corev
      IDATAN(2)=nucev
      IDATAN(3)=nuccv

      CNAMEN='NEUTVER'

      NWMASTN=3
C      NSEGMN=corev
      NSEGMN=1
      LSEGMN=3

      CALL KZBCR1(CNAMEN,NWMASTN,'3I',NSEGMN,LSEGMN,'3I',IERR)      
      IF (IERR.NE.0) GOTO 1000
      CALL KZREP1(CNAMEN,0,NWMASTN,IDATAN(1))


C--- Flux and geometry
      IDATAF(1)=SK_GEOMETRY
      IDATAF(2)=0
      RDATAF(3)=0
      IDATAF(4)=0
      RDATAF(5)=0
      IDATAF(6)=ITAUFLG

      CNAMEF='NEUTFLUX'

      NWMASTF=6
C      NSEGMF=GEOMVER
      NSEGMF=1
      LSEGMF=6

      CALL KZBCR1(CNAMEF,NWMASTF,'6I',NSEGMF,LSEGMF,'6I',IERR)
      IF (IERR.NE.0) GOTO 1000
      CALL KZREP1(CNAMEF,0,NWMASTF,IDATAF(1))


C--- Configuration of nuclear efffects, etc...
      IDATAC(1)=NEFRMFLG
      IDATAC(2)=NEPAUFLG
      IDATAC(3)=NENEFO16
      IDATAC(4)=NEMODFLG
      RDATAC(5)=PFSURF
      RDATAC(6)=PFMAX
      RDATAC(7)=VNUINI
      RDATAC(8)=VNUFIN
      IDATAC(9)=IFORMLEN
      IDATAC(10)=NUCRESCAT
      RDATAC(11)=XNUCFACT

      CNAMEC='NEUTCONFIG'

      NWMASTC=4
C      NSEGMC=NEFRMFLG
      NSEGMC=1
      LSEGMC=4

      CALL KZBCR1(CNAMEC,NWMASTC,'4I',NSEGMC,LSEGMC,'4I',IERR)      
      IF (IERR.NE.0) GOTO 1000
      CALL KZREP1(CNAMEC,0,NWMASTC,IDATAC(1))


C--- Reaction model
      IDATAM(1)=MODELDIS
      IDATAM(2)=MODELCOH

      CNAMEM='NEUTMODEL'

      NWMASTM=2
C      NSEGMM=MODELDIS
      NSEGMM=1
      LSEGMM=7

      CALL KZBCR1(CNAMEM,NWMASTM,'2I',NSEGMM,LSEGMM,'5I2F',IERR)      
      IF (IERR.NE.0) GOTO 1000
      CALL KZREP1(CNAMEM,0,NWMASTM,IDATAM(1))

      IDATAM(1)=MDLQE
      IDATAM(2)=MDLSPI
      IDATAM(3)=MDLCOH
      IDATAM(4)=MDLDIS
      IDATAM(5)=MDLQEAF
      RDATAM(6)=XMAQE
      RDATAM(7)=XMASPI
      CALL KZREP1(CNAMEM,1,LSEGMM,IDATAM(1))      

#ifdef VERBOSE
      print*, 'NEUTVER-> ',IDATAN
      print*, 'NEUTFLUX-> ',
     &     IDATAF(1),IDATAF(2),RDATAF(3),IDATAF(4),
     &     RDATAF(5),IDATAF(6)
      print*, 'NEUTCONFIG-> ',
     &     IDATAC(1),IDATAC(2),IDATAC(3),IDATAC(4),
     &     RDATAC(5),RDATAC(6),RDATAC(7),RDATAC(8),IDATAC(9),
     &     IDATAC(10),RDATAC(11)
      print*, 'NEUTMODEL-> ',IDATAM
#endif

      RETURN

C     
C     ++ ERROR MESSAGE
C     
 1000 WRITE(6,900)
  900 FORMAT(' *** ERROR *** ( BANK CREATION IN NEMKMODELBK )')
      RETURN
      END

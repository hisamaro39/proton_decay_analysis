***********************************************************************
*     ------------------------
      SUBROUTINE NEMKNEBK(POS)
*     ------------------------
*     
*     ( purpose )
*       Create NEUT bank
*
*     ( input ) 
*       POS : INTERACTION POINT
*       COMMON /NEWORK/
*       COMMON /NESOLACT/      
*
*     ( output ) 
*       NEUT BANK
*
*     ( creation date and author )
*       1995.03.15 ; First version by K.Kaneyuki
*
***********************************************************************
#include "nework.h"
#include "mcgenpar.h"      
      REAL    POS(3)
      INTEGER IDATA(500)
      REAL    RDATA(500)
      CHARACTER*32 CNAME
      EQUIVALENCE (IDATA(1),RDATA(1))
C
C     -- MAKE EMB
C
C     LEM+1 : # OF PARTICLE
C     LEM+2 : MODE OF INTERACTION
C     LEM+3 : ISOL   ( FLAG OF SOLACT )
C     LEM+4 : SOLACT ( SOLAR ACTIVITY )
C
      IDATA(1)=NUMNE
      IDATA(2)=MODENE
      IDATA(3)=0
      RDATA(4)=0
      RDATA(5)=POS(1)
      RDATA(6)=POS(2)
      RDATA(7)=POS(3)

      CNAME='NEUT'

      NWMAST=7
      NSEGM=NUMNE
      LSEGM=7

      CALL KZBCR1(CNAME,NWMAST,'3I 4F',NSEGM,LSEGM,'1I 3F 3I',IERR)      
      IF (IERR.NE.0) GOTO 1000
      CALL KZREP1(CNAME,0,NWMAST,IDATA(1))      

      DO 10 I=1,NUMNE
C     
C     -- SET VECTOR INFORMATION OF M.C. ON I-TH PARTICLE
C     
C        IND+ 1 : PARTICLE CODE
C           + 2 : PX ( MEV )
C           + 3 : PY ( MEV )
C           + 4 : PZ ( MEV )
C           + 5 : PARENT PARTICLE
C           + 6 : FLAG OF FINAL STATE
C           + 7 : CHERENKOV FLAG
C
         IDATA( 1)=IPNE(I)
         RDATA( 2)=PNE(1,I)
         RDATA( 3)=PNE(2,I)
         RDATA( 4)=PNE(3,I)
         IDATA( 5)=IORGNE(I)
         IDATA( 6)=IFLGNE(I)
         IDATA( 7)=ICRNNE(I)
C     
C     -- MAKE I-TH SEGMENT
C     
         ISEG=I
         CALL KZREP1(CNAME,ISEG,7,IDATA)         

 10   CONTINUE

      RETURN
C     
C     ++ ERROR MESSAGE
C     
 1000 WRITE(6,900)
  900 FORMAT(' *** ERROR *** ( BANK CREATION IN NEMKNEBK )')
      RETURN
      END

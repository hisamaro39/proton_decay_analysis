***********************************************************************
*     ------------------------
      SUBROUTINE NERDCRSBK
*     ------------------------
*     
*     ( purpose )
*       Read NEUTCRS bank
*
*     ( input ) 
*       none
*
*     ( output ) 
*       COMMON /NEUTCRS/      
*
*     ( creation date and author )
*       1995.03.16 ; First version by K.Kaneyuki
*       2011.05.13 ; Modified for NEUTCRS bank by P.de Perio
*
***********************************************************************
#include "neutcrs.h"

      CHARACTER*32 CNAME
      INTEGER IDATA(20)
      REAL    RDATA(20)
      EQUIVALENCE (RDATA(1),IDATA(1))

      CNAME='NEUTCRS'
      CALL KZBLOC(CNAME,IPBANK)
      IF (IPBANK.EQ.0) GOTO 1000
C
C     -- READ EMB FOR VERTEX BANK
C
C     LEM+1 : # OF VERTICES
C
      ISEG=0
      CALL KZGET1(CNAME,ISEG,NDATA,RDATA)
      IF (NDATA.EQ.0) GOTO 2000
      CRSX   = RDATA(1)
      CRSY   = RDATA(2)
      CRSZ   = RDATA(3)
      CRSPHI = RDATA(4)
      
      RETURN
C     
C     ++ ERROR MESSAGE
C     
 1000 WRITE(6,900)
 900  FORMAT(' *** ERROR *** ( NEUTCRS DOES NOT EXIST IN NERDCRSBK )')
      RETURN
 2000 WRITE(6,910) ISEG
 910  FORMAT(' *** ERROR *** ( NEUTCRS BANK : ',I2,
     $     'DOES NOT EXIST IN NERDCRSBK )')
      RETURN
      END

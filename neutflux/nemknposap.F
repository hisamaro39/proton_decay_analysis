***********************************************************************
*     -----------------
      SUBROUTINE nemknposap
*     -----------------
*
*     ( purpose )
*       Create NUCPOS bank
*
*     ( input ) 
*       none
*
*     ( output ) 
*       none
*
*     ( creation date and author )
*       1992.6.5 ; First version by K.Kaneyuki
*       1994.9.19; mod. for zbs V1.5
*       2006.6.12; Store all the positions
*
***********************************************************************
      IMPLICIT NONE

      INTEGER IDATA(20)
      REAL    RDATA(20)
      EQUIVALENCE(IDATA(1),RDATA(1))
      CHARACTER*32 CNAME

      integer*4 ISEG,NWMAST,LSEGM,NSEGM,IERR
      integer*4 IPBANK,NDATA

#include "vcwork.h"
#include "nework.h"
#include "posinnuc.h"

      CNAME='NUCPOS'

      CALL KZBLOC(CNAME,IPBANK)
      IF (IPBANK.NE.0) GOTO 5000

C
C     -- READ EMB
C
C     LEM+1 : # OF VERTEX
C

      ISEG = 0

      CALL KZGET1(CNAME,ISEG,NDATA,IDATA)
      IF (NDATA.NE.1) GOTO 2000
      NSEGM = IDATA(1)

      DO 200 ISEG=1,NSEGM
C     
C     -- GET LOCATION in the nucleus INFORMATION 
C     
C        IND+ 1 : X
C           + 2 : Y
C           + 3 : Z
C     
C     -- MAKE I-TH SEGMENT
C     
C         write(*,'(A,F12.5,A,F12.5,A,F12.5,A)')
C     $        'Position:(',posnuc(1,ISEG),
C     $                 ',',posnuc(2,ISEG),
C     $                 ',',posnuc(3,ISEG),')'
C
         CALL KZGET1(CNAME,ISEG,NDATA,IDATA)

         IF (NDATA.ne.3) goto 3000
         
         RDATA( 1)=POSNUC(1,ISEG)
         RDATA( 2)=POSNUC(2,ISEG)
         RDATA( 3)=POSNUC(3,ISEG)

 200  continue

      RETURN
C     
C     ++ ERROR MESSAGE
C     
 2000 WRITE(6,2900)
 2900 FORMAT('nerdnpos: ERR :# of position was not properly stored.')
      goto 5000

 3000 WRITE(6,3900) 'nerdnpos: ERR :# of data is not 3.',
     $     NDATA,' at ',ISEG
 3900 FORMAT(6(A30,I4,A5,I4))

 5000 RETURN
      END

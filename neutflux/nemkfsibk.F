***********************************************************************
*     ------------------------
      SUBROUTINE NEMKFSIBK
*     ------------------------
*     
*     ( purpose )
*       Create NEUT FSIHIST bank
*
*     ( input ) 
*       COMMON /FSIHIST/
*
*     ( output ) 
*       NEUT BANK
*
*     ( creation date and author )
*       1995.03.16 ; First version by K.Kaneyuki
*       2010.10.24 ; Modified for FSIHIST bank by P.de Perio
*
***********************************************************************
#include "fsihist.h"

      INTEGER IDATA(500)
      REAL    RDATA(500)
      CHARACTER*32 CNAME
      EQUIVALENCE (RDATA(1),IDATA(1))
C
C     -- MAKE EMB FOR VERTEX BANK
C
C     LEM+1 : # OF VERTICES

      IDATA(1)=NVERT

      CNAME='FSIHISTVERT'

      NWMAST=1
      NSEGM=NVERT
      LSEGM=4

      CALL KZBCR1(CNAME,NWMAST,'1I',NSEGM,LSEGM,'3F 1I',IERR)      
      IF (IERR.NE.0) GOTO 1000
      CALL KZREP1(CNAME,0,NWMAST,IDATA(1))      

      DO 10 I=1,NVERT
C     
C     -- SET VECTOR INFORMATION OF M.C. ON I-TH VERTEX
C     
C        IND+ 1 : POS X (fm)
C           + 2 : POS Y (fm)
C           + 3 : POS Z (fm)
C           + 4 : INTERACTION FLAG OF VERTEX (see fsihist.h)
C
         RDATA( 1)=POSVERT(1,I)
         RDATA( 2)=POSVERT(2,I)
         RDATA( 3)=POSVERT(3,I)
         IDATA( 4)=IFLGVERT(I)
C     
C     -- MAKE I-TH SEGMENT
C     
         ISEG=I
         CALL KZREP1(CNAME,ISEG,LSEGM,RDATA)         

 10   CONTINUE


C
C     -- MAKE EMB FOR INTERMEDIATE PARTICLE BANK
C
C     LEM+1 : # OF INTERMEDIATE PARTICLES

      IDATA(1)=NVCVERT
      RDATA(2)=FSIPROB

      CNAME='FSIHISTPART'

      NWMAST=2
      NSEGM=NVCVERT
      LSEGM=8

      CALL KZBCR1(CNAME,NWMAST,'1I 1F',NSEGM,LSEGM,'5F 3I',IERR)      
      IF (IERR.NE.0) GOTO 1000
      CALL KZREP1(CNAME,0,NWMAST,IDATA(1))      

      DO 20 I=1,NVCVERT
C     
C     -- SET VECTOR INFORMATION OF M.C. ON I-TH INTERMEDIATE PARTICLE
C     
C        IND+ 1 : DIR X
C           + 2 : DIR Y
C           + 3 : DIR Z
C           + 4 : ABSOLUTE MOM. in lab frame (MeV/c)
C           + 5 : ABSOLUTE MOM. in nucleon rest frame (MeV/c)
C           + 6 : PDG Particle code
C           + 7 : Index of initial vertex
C           + 8 : Index of final vertex
C
         RDATA( 1)=DIRVERT(1,I)
         RDATA( 2)=DIRVERT(2,I)
         RDATA( 3)=DIRVERT(3,I)
         RDATA( 4)=ABSPVERT(I)
         RDATA( 5)=ABSTPVERT(I)
         IDATA( 6)=IPVERT(I)
         IDATA( 7)=IVERTI(I)
         IDATA( 8)=IVERTF(I)
C     
C     -- MAKE I-TH SEGMENT
C     
         ISEG=I
         CALL KZREP1(CNAME,ISEG,LSEGM,RDATA)         

 20   CONTINUE


      RETURN
C     
C     ++ ERROR MESSAGE
C     
 1000 WRITE(6,900)
  900 FORMAT(' *** ERROR *** ( BANK CREATION IN NEMKFSIBK )')
      RETURN
      END

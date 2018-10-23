***********************************************************************
*     ------------------------
      SUBROUTINE NERDFSIBK
*     ------------------------
*     
*     ( purpose )
*       Read NEUT FSIHIST bank
*
*     ( input ) 
*       none
*
*     ( output ) 
*       COMMON /FSIHIST/      
*
*     ( creation date and author )
*       1995.03.16 ; First version by K.Kaneyuki
*       2010.10.24 ; Modified for FSIHIST bank by P.de Perio
*
***********************************************************************
#include "fsihist.h"

      CHARACTER*32 CNAME
      INTEGER IDATA(20)
      REAL    RDATA(20)
      EQUIVALENCE (RDATA(1),IDATA(1))

      CNAME='FSIHISTVERT'
      CALL KZBLOC(CNAME,IPBANK)
      IF (IPBANK.EQ.0) GOTO 1000
C
C     -- READ EMB FOR VERTEX BANK
C
C     LEM+1 : # OF VERTICES
C
      ISEG=0
      CALL KZGET1(CNAME,ISEG,NDATA,IDATA)
      IF (NDATA.EQ.0) GOTO 2000
      NVERT =IDATA(1)

      DO 10 I=1,NVERT
C     
C     -- READ VECTOR INFORMATION OF M.C. ON I-TH VERTEX
C
C        IND+ 1 : POS X (fm)
C           + 2 : POS Y (fm)
C           + 3 : POS Z (fm)
C           + 4 : INTERACTION FLAG OF VERTEX (see fsihist.h)
C
         ISEG=I
         CALL KZGET1(CNAME,ISEG,NDATA,RDATA)

         POSVERT(1,I) =RDATA(1)
         POSVERT(2,I) =RDATA(2)
         POSVERT(3,I) =RDATA(3)
         IFLGVERT(I)  =IDATA(4)

 10   CONTINUE



      CNAME='FSIHISTPART'
      CALL KZBLOC(CNAME,IPBANK)
      IF (IPBANK.EQ.0) GOTO 1000
C
C     -- READ EMB FOR INTERMEDIATE PARTICLE BANK
C
C     LEM+1 : # OF INTERMEDIATE PARTICLES
C
      ISEG=0
      CALL KZGET1(CNAME,ISEG,NDATA,IDATA)
      IF (NDATA.EQ.0) GOTO 2000
      NVCVERT =IDATA(1)
      FSIPROB =RDATA(2)

      DO 20 I=1,NVCVERT
C     
C     -- READ VECTOR INFORMATION OF M.C. ON I-TH INTERMEDIATE PARTICLE
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
         ISEG=I
         CALL KZGET1(CNAME,ISEG,NDATA,RDATA)

         DIRVERT(1,I) =  RDATA( 1)
         DIRVERT(2,I) =  RDATA( 2)
         DIRVERT(3,I) =  RDATA( 3)
         ABSPVERT(I)  =  RDATA( 4)
         ABSTPVERT(I) =  RDATA( 5)
         IPVERT(I)    =  IDATA( 6)
         IVERTI(I)    =  IDATA( 7)
         IVERTF(I)    =  IDATA( 8)

 20   CONTINUE



      RETURN
C     
C     ++ ERROR MESSAGE
C     
 1000 WRITE(6,900)
 900  FORMAT(' *** ERROR *** ( NEUT DOES NOT EXIST IN NERDFSIBK )')
      RETURN
 2000 WRITE(6,910) ISEG
 910  FORMAT(' *** ERROR *** ( NEUT BANK : ',I2,
     $     'DOES NOT EXIST IN NERDFSIBK )')
      RETURN
      END

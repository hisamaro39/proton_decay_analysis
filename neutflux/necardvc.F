************************************************************************
*     -------------------
      SUBROUTINE NECARDVC
*     -------------------
*
*     (Purpose)
*       Read VECT card and setup common blocks.
*
*     (Input)
*       card file
*
*     (Output)
*
*     (Creation Date and Author)
*       1995.1.30 ; Koshio ( SGCARD )
*       1995.11.17 ; K.Kaneyuki
*
************************************************************************
#include "skhead.h"
#include "dsbnkf.h"
#include "necardvc.h"

      CHARACTER*40 TMPFIL, COMMAND
      INTEGER*4 GEOMVER

C --  Initialize

      NERDFLG = 0
      NEVENTCD = 0
      NEPARTCD = 0
      NEVTXFLG = 0
      CALL VZERO(3,POSCRDNE)
      DWALLNE = 0.0
      NEDIRFLG = 0
      CALL VZERO(3,DIRCRDNE)
      NEMOMFLG = 0
      CALL VZERO(2,PCRDNE)
      SK_GEOMETRY = 0
      
C --- get card file name (1st argument of the program)

      WRITE(6,*)
      IF (IARGC() .GE. 1) THEN
         CALL GETARG(1, FCARD)
      ELSE
         FCARD = 'supersim.card'
      END IF
      WRITE(6,*) '======= CONFIGURATION FOR ATMPD VECT ======='
      WRITE(6,*) 'NEUT(VECT) card file you read is ',FCARD

C -- extract NEUT card

      CALL GETUNI(TMPFIL)
      OPEN(UNIT=90, FILE=FCARD, STATUS='OLD')
      OPEN(UNIT=91, FILE=TMPFIL,STATUS='UNKNOWN')
      CALL EXCARD('VECT', 90, 91, IRET)
      CLOSE(UNIT=91)
      CLOSE(UNIT=90)

      OPEN(UNIT=91, FILE=TMPFIL,STATUS='OLD')
      CALL FFINIT(0)
      CALL FFSET('LINP', 91)
      CALL FFSET('SIZE', 10)
      CALL FFKEY('FILE',NERDFLG, 1,'INTEGER')
      CALL FFKEY('NEVT',NEVENTCD,1,'INTEGER')
      CALL FFKEY('PART',NEPARTCD,1,'INTEGER')
      CALL FFKEY('MVTX',NEVTXFLG,1,'INTEGER')
      CALL FFKEY('VTX' ,POSCRDNE,3,'REAL')
      CALL FFKEY('WALL',DWALLNE, 1,'REAL')
      CALL FFKEY('MDIR',NEDIRFLG,1,'INTEGER')
      CALL FFKEY('DIR' ,DIRCRDNE,3,'REAL')
      CALL FFKEY('MMOM',NEMOMFLG,1,'INTEGER')
      CALL FFKEY('MOM' ,PCRDNE  ,2,'REAL')
      CALL FFKEY('GEOMVER' ,GEOMVER  ,1,'INTEGER')
      CALL FFGO
      CLOSE(UNIT=91)
      COMMAND = '/bin/rm '//TMPFIL
      CALL SYSTEM(COMMAND)

      WRITE(6,*) 'NECARDVC : INPUTTED PARAMETER '
      WRITE(6,*) 'NECARDVC : FILE = ',NERDFLG
      WRITE(6,*) 'NECARDVC : NEVT = ',NEVENTCD
      IF (NERDFLG.EQ.2) THEN
         WRITE(6,*) 'NECARDVC : PART = ',NEPARTCD
         WRITE(6,*) 'NECARDVC : MVTX = ',NEVTXFLG
         IF (NEVTXFLG.EQ.1) WRITE(6,*) 'NECARDVC : VTX = ',POSCRDNE
         IF (NEVTXFLG.EQ.2) WRITE(6,*) 'NECARDVC : WALL = ',DWALLNE
         WRITE(6,*) 'NECARDVC : MDIR = ',NEDIRFLG
         IF (NEDIRFLG.EQ.1) WRITE(6,*) 'NECARDVC : DIR = ',DIRCRDNE
         WRITE(6,*) 'NECARDVC : MMOM = ',NEMOMFLG
         WRITE(6,*) 'NECARDVC : MOM  = ',PCRDNE
      ENDIF
      WRITE(6,*) 'NECARDVC : GEOMVER(GEOMETRY VERSION) = ',GEOMVER
      if (GEOMVER.le.0) then
         write(6,*) 'GEOMVER must be set.'
         STOP
      endif
      SK_GEOMETRY = GEOMVER

      RETURN
      END



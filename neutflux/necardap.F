************************************************************************
*     -----------------
      SUBROUTINE NECARDAP
*     -----------------
*
*     (Purpose)
*       Read NEUT card and setup common blocks.
*
*     (Input)
*       card file
*
*     (Output)
*
*     (Creation Date and Author)
*       1995.1.30 ; Koshio ( SGCARD )
*       1995.11.17 ; K.Kaneyuki
*       1999.10.11 ; D.Casper - allow user to specify #years in card file
*       2003.03.30 ; add neutparams.h(PFSURF,PFMAX,VNUINI,VNUFIN)
*       2007.04.24 ; G.Mitsuka necard->necardap(for ATMPD)
*                    add NECOHEPI
*       2007.12.05 ; G.Mitsuka added NEWATERVOL (default is 0)
************************************************************************

#include "dsbnkf.h"
#include "necardap.h"
#include "nesolact.h"      
#include "neutparams.h"

      CHARACTER*40 TMPFIL, COMMAND

C --  Initialize

      SOLACT = 0.
      NEFLXFLG = 0
      NEYEARS = 1.0
      ITAUFLG = 0
      NEWATERVOL = 0
      
C --- get card file name (1st argument of the program)

      WRITE(6,*) 
      IF (IARGC() .GE. 1) THEN
         CALL GETARG(1, FCARD)
      ELSE
         FCARD = 'neut.card'
      END IF
      WRITE(6,*) '======= CONFIGURATION FOR ATMPD ======='
      WRITE(6,*) 'NEUT card file you read is ',FCARD

C -- extract NEUT card

      CALL GETUNI(TMPFIL)
      OPEN(UNIT=90, FILE=FCARD, STATUS='OLD')
      OPEN(UNIT=91, FILE=TMPFIL,STATUS='UNKNOWN')
      CALL EXCARD('NEUT', 90, 91, IRET)
      CLOSE(UNIT=91)
      CLOSE(UNIT=90)

      OPEN(UNIT=91, FILE=TMPFIL,STATUS='OLD')
      CALL FFINIT(0)
      CALL FFSET('LINP', 91)
      CALL FFSET('SIZE', 9)

      CALL FFKEY('ISOL',ISOL    , 1,'INTE')
      CALL FFKEY('SOL' ,SOLACT  , 1,'REAL')
      CALL FFKEY('FLX' ,NEFLXFLG, 1,'INTE')
      CALL FFKEY('YEAR',NEYEARS,1,'REAL')
      CALL FFKEY('TAU',ITAUFLG,1,'INTEGER')
      CALL FFKEY('WATERVOL',NEWATERVOL,1,'INTEGER')

      CALL FFGO
      CLOSE(UNIT=91)
      COMMAND = '/bin/rm '//TMPFIL
      CALL SYSTEM(COMMAND)

      WRITE(6,*) 'NECARDAP : INPUTTED PARAMETER '
      WRITE(6,*) 'NECARDAP : ISOL = ',ISOL
      IF (ISOL.EQ.0) THEN
         SOLACT = 0.5
      ELSE
         WRITE(6,*) 'NECARDAP : SOL  = ',SOLACT
      ENDIF
      WRITE(6,*) 'NECARDAP : NEFLXFLG = ',NEFLXFLG
      WRITE(6,*) 'NECARDAP : Generating ',NEYEARS,' years of data'
      WRITE(6,*) 'NECARDAP : RANDOM= ',NEYEARS,' years of data'
      WRITE(6,*) 'NECARDAP : TAU = ',ITAUFLG
      WRITE(6,*) 'NECARDAP : WATERVOLUME = ',NEWATERVOL

      RETURN
      END



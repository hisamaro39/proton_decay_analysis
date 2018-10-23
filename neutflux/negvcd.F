************************************************************************
*     ------------------------------
      SUBROUTINE NEGVCD(ID,POS,PMOM)
*     ------------------------------
*
*     (Purpose)
*        Copy vector information for neut
C
C     (INPUT)
C        COMMON /NECARD2/
C
C     (OUTPUT)
*        ID   : Particle ID
C        POS  : Position
C        PMOM : Momentum
*        
C     (Creation Date and Author)
C       199?.??.?? ; ????
*
************************************************************************
#include "necardvc.h"

      REAL POS(3),PMOM(3)
C
C --  SET ID,POSITION,MOMMENTUM
C
      ID=NEPARTCD

      IF(NEVTXFLG.EQ.1) THEN
         POS(1)=POSCRDNE(1)
         POS(2)=POSCRDNE(2)
         POS(3)=POSCRDNE(3)
      ELSE
         IF(NEVTXFLG.EQ.2) CALL RNPOS(POS,DWALLNE)
      ENDIF

      IF(NEDIRFLG.EQ.1) THEN
         PMOM(1)=DIRCRDNE(1)
         PMOM(2)=DIRCRDNE(2)
         PMOM(3)=DIRCRDNE(3)
      ELSE
         IF(NEDIRFLG.EQ.2) CALL RNDIR(PMOM)
      ENDIF

      IF(NEMOMFLG.EQ.1) THEN
         PP = PCRDNE(1)
      ELSE
         IF(NEMOMFLG.EQ.2) PP = PCRDNE(1)+
     &                          RLU(DUMMY)*(PCRDNE(2)-PCRDNE(1))
      ENDIF
      PMOM(1)=PMOM(1) * PP
      PMOM(2)=PMOM(2) * PP
      PMOM(3)=PMOM(3) * PP

      RETURN
      END

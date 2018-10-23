***********************************************************************
*     ------------------
      SUBROUTINE NERDLEE
*     ------------------
*
*     ( Purpose )
*       READ FLUX CALCULATED BY LEE, AND FILL INTO COMMON / NELEEF /
*
*     ( Creation Date and Auther )
*       1991.07.01 ; T.KAJITA  
*       1996.01.23 ; Y.KANAYA
*                    FOR SK
*
******************************************************************
#include "skcnst.h"
      PARAMETER(LUN=71)
      REAL    ENU(136),DENU(135)
      REAL    FLUX(4,136,15),FLXTOT(4,136)
      CHARACTER*1  BUFF1
      COMMON / NELEEF / ENELOW(135),
     &       FEMIN(135), FEMAX(135), FEBMIN(135), FEBMAX(135),
     &       FMMIN(135), FMMAX(135), FMBMIN(135), FMBMAX(135)
C
      CALL VZERO(FLXTOT,544)
C
      CALL SKOPENF(LUN,1,'F',IERR)
      IF(IERR.NE.0) STOP      
C
      DO 10 I=1,4
        DO 20 K=1,15
          READ(LUN,'(A)') BUFF1
          READ(LUN,'(A)') BUFF1
          DO 30 J=1,27
            READ(LUN,*) (FLUX(I,J*5-5+L,K),L=1,5)
30        CONTINUE
          READ(LUN,*) FLUX(I,136,K)
20      CONTINUE
10    CONTINUE
C
      CALL SKCLOSEF(NFILE)
C
      DO 1900 I=1,4
      DO 1901 J=1,136
      DO 1902 K=1,15
      FLXTOT(I,J)=FLXTOT(I,J)+FLUX(I,J,K)
 1902 CONTINUE
 1901 CONTINUE
 1900 CONTINUE
C
      DO 40 I=1,50
        ENU(I)=(I-1)*5.+2.5
40    CONTINUE

      DO 50 I=1,25
        ENU(I+50)=(I-1)*10.+255.
50    CONTINUE

      DO 60 I=1,25
        ENU(I+75)=(I-1)*20.+510.
60    CONTINUE

      DO 70 I=1,20
        ENU(I+100)=(I-1)*50.+1025.
70    CONTINUE

      DO 80 I=1,15
        ENU(I+120)=(I-1)*200+2100.
80    CONTINUE

      ENU(136)=5100.

      DO 90 I=1,50
        DENU(I)=5.
90    CONTINUE

      DO 91 I=51,75
        DENU(I)=10.
91    CONTINUE

      DO 92 I=76,100
        DENU(I)=20.
92    CONTINUE

      DO 93 I=101,120
        DENU(I)=50.
93    CONTINUE

      DO 94 I=121,135
        DENU(I)=200.
94    CONTINUE
C
C ++ NORMARIZE TO /M**2/SEC/SR/GEV
      DO 2900 I=1,4
*      DO 2901 J=1,136
      DO 2901 J=1,135
      FLXTOT(I,J)=FLXTOT(I,J)*3.981E-3/(4.*PI)/(DENU(J)/1000.)
 2901 CONTINUE
 2900 CONTINUE
C
C
C COPY TO COMMOM /NELEEF/
C
      DO 2000 J=1,135
      ENELOW(J)= ENU(J) /1000.
      FEMIN(J)=FLXTOT(1,J)
      FEBMIN(J)=FLXTOT(2,J)
      FMMIN(J)=FLXTOT(3,J)
      FMBMIN(J)=FLXTOT(4,J)
      FEMAX(J)=FLXTOT(1,J)
      FEBMAX(J)=FLXTOT(2,J)
      FMMAX(J)=FLXTOT(3,J)
      FMBMAX(J)=FLXTOT(4,J)
 2000 CONTINUE
C
      RETURN
      END

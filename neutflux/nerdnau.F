*********************************************************************
*     ------------------
      SUBROUTINE NERDNAU
*     ------------------
*
*     ( Purpose )
*       READ FLUX CALCULATED BY NUAMOV, AND FILL INTO COMMON /NENAUMF/
*
*     ( Creation Date and Auther )
*       1992.09.03 ; T.KAJITA  
*       1996.01.23 ; Y.KANAYA
*                    FOR SK
*
*********************************************************************
      PARAMETER(LUN=72)
      DIMENSION FLINT(4),FLDIF(4)
      COMMON /NENAUMF / ENEU(33),
     &       FEMIN(33), FEMAX(33), FEBMIN(33), FEBMAX(33),
     &       FMMIN(33), FMMAX(33), FMBMIN(33), FMBMAX(33)
C
      CALL SKOPENF(LUN,1,'F',IERR)
      IF(IERR.NE.0) STOP      
C
      READ(LUN,'(A)') BUFF1
      READ(LUN,'(A)') BUFF1
      READ(LUN,'(A)') BUFF1
C
      DO 10 IE=1,33
*      READ(LUN,500,END=100) EN,(FLINT(I),I=1,4),(FLDIF(J),J=1,4),
      READ(LUN,*) EN,(FLINT(I),I=1,4),(FLDIF(J),J=1,4),
     & FNENM, RNMBNM, RNEBNE
* 500  FORMAT(E10.4,8F10.3,3F6.3)
CC    WRITE(6,510) IE,PCOSL(ICOS),EN,(FLINT(I),I=1,4),
CC   &             (FLDIF(J),J=1,4), FNENM, RNMBNM, RNEBNE
CC510 FORMAT(I3,F5.2,2X,E10.4,4F8.2,4F9.2,3F6.3)
C
      CALL SKCLOSEF(NFILE)
C
C      ENELOW(IE)=EN
      ENEU(IE)=EN
      FEMIN(IE)=FLDIF(3)
      FEBMIN(IE)=FLDIF(4)
      FMMIN(IE)=FLDIF(1)
      FMBMIN(IE)=FLDIF(2)
      FEMAX(IE)=FLDIF(3)
      FEBMAX(IE)=FLDIF(4)
      FMMAX(IE)=FLDIF(1)
      FMBMAX(IE)=FLDIF(2)
C
  10  CONTINUE
 100  CONTINUE
C
      RETURN
      END

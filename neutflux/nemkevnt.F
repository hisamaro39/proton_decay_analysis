**********************************************************************
*     ----------------------------
#ifdef UPMU_RCK
      SUBROUTINE NEMKEVNT(IPAR,UPMU_FLAG)
#else
      SUBROUTINE NEMKEVNT(IPAR,POS)
#endif
*     ----------------------------
*     
*     ( purpose )
*       VECTOR GENERATION ( NEUTRINO INTERACTION )
*
*     ( input ) 
*       IPAR    : PARICLE CODE
*       POS(3)  : Vertex position
*
*     ( output ) 
*       COMMON NEWORK
*
*     ( creation date and author )
*       1995.01.12 ; First version by Y.Hayato
*       2007.11.10 ; ADD UPMU MODE by T.Tanaka      
*
**********************************************************************
#ifdef UPMU_RCK
      IMPLICIT NONE
#include "upmu.h"
#include "geotnk.h"
      REAL*4 THETA,PHI
*---- *** get muon info.
      INTEGER*4 NW,NSEG,ISEG,I,IDATA(20),IP,
     $     ISTA_NU_UPMU, ISTA_MU_UPMU
      REAL*4 RDATA(20),INF_POS(3),MUPM,INF_PMOM(3),PM,
     $      MUDIR(3),MUPOS(3),NUDIR(3),DIR(3),ENEUGEV
      REAL*4 V3ABS,V3COS
      INTEGER*4 UPMU_FLAG
      EXTERNAL V3ABS,V3COS
      EQUIVALENCE(IDATA(1),RDATA(1))
#endif

      INTEGER*4 IPAR
      REAL*4 POS(3)
      
      REAL PMOM(3),DIRNEU(3)

      EXTERNAL RNNEUE
      REAL*4 RNNEUE

      REAL*4 ENEU

#ifdef UPMU_RCK
      UPMU_FLAG=0
*----- PID check
      IF(ABS(IPAR).NE.14) THEN
         WRITE(6,*) "PARTICLE BAD!!!",IPAR
         STOP
      ENDIF
*----- GENERATE
 5    ENEUGEV=RNNEUE(IPAR)      !---- energy(GeV)
      CALL NERNNEUT(ENEUGEV,DIRNEU,IPAR) !--- dir.
      CALL DIRTRANS(DIRNEU,THETA,PHI)    !--- trans dir.
      CALL RNPOS_RCK(POS,RROCK,THETA,PHI)!--- rock position
      ENEU=ENEUGEV*1.0E3        !---- energy(MeV)
      DO 10 I=1,3
         PMOM(I)=ENEU*DIRNEU(I)
 10         CONTINUE
      CALL CHK_NU_UPMU(POS,ENEU,DIRNEU,ISTA_NU_UPMU) !-- chk. nU
      IF(ISTA_NU_UPMU.EQ.1) THEN
          WRITE(6,*) IPAR, POS, PMOM
          CALL NEVECTAP(IPAR,POS,PMOM,IERR) !- make mu vector
          if (ierr.ne.0) then
             goto 5
          endif
*---  *** get muon info.
*     *** get input particle number
         call kzget1('MCVECT',0,nw,idata) !-- get vec. info.
         nseg = idata(1)                  !-- # of particle

*     *** vertex bank
         call kzget1('MCVERTEX',1,nw,idata) !-- get nu info.
         inf_pos(1)=rdata(1)
         inf_pos(2)=rdata(2)
         inf_pos(3)=rdata(3)

*     *** loop for input particle
         do 300 iseg= 1,nseg
            call kzget1('MCVECT',iseg,nw,idata)
*     *** input data
            ip = idata(1)
            inf_pmom(1) = rdata(3)
            inf_pmom(2) = rdata(4)
            inf_pmom(3) = rdata(5)
            pm = v3abs(inf_pmom)
            if(pm.eq.0) goto 300
            dir(1)=inf_pmom(1)/pm
            dir(2)=inf_pmom(2)/pm
            dir(3)=inf_pmom(3)/pm
            if(abs(ip).eq.14) then
              CALL v3copy(dir,nudir)
            elseif(abs(ip).eq.13) then
               CALL v3copy(dir,mudir)
               CALL v3copy(inf_pos,mupos)
               mupm=pm
            endif
 300     continue
         if(mudir(1).ne.0.and.mudir(2).ne.0.
     &        and.mudir(3).ne.0) then
              
*---  *** store check
            CALL CHK_MU_UPMU(MUPOS,MUPM,MUDIR,ISTA_MU_UPMU)
            IF(ISTA_MU_UPMU.EQ.1) UPMU_FLAG=1
         ENDIF
      ENDIF
#else


 5    IERR = 0
      ENEU=RNNEUE(IPAR)
      CALL NERNNEUT(ENEU,DIRNEU,IPAR)

      ENEU=ENEU*1.0E3
      DO 10 I=1,3
         PMOM(I)=ENEU*DIRNEU(I)
 10   CONTINUE

      CALL NEVECTAP(IPAR,POS,PMOM,IERR)
      if (ierr.ne.0) then
         goto 5
      endif
#endif
      RETURN
      END
      

C     
C     TEST PROGRAM TO READ VECTOR
C
      IMPLICIT NONE


#include "nework.h"
#include "vcwork.h"
#include "vcvrtx.h"

C#include "skhead.h"
C#include "skparm.h"
C#include "sktq.h"

      INTEGER*4 LUNI
      PARAMETER(LUNI=10)

      REAL POS(3)

      INTEGER*4  LRECL
      PARAMETER (LRECL=1024)

      INTEGER*4  ISIZE
      PARAMETER (ISIZE=5000000)

      REAL*4     H(ISIZE)
      COMMON/PAWC/H

      CHARACTER*128 FNAME_NTPL

      INTEGER*4  IERR,I,N1,N2,ICYCLE,LBANK

C     Neut PARAMETER COMMON
      INTEGER*4   IEVENT,MODEN,NPARNEUT,IPNEUT(100)
      REAL*4      ABSPNEUT(100),PNEUT(3,100)
      COMMON/COM1/IEVENT,MODEN,NPARNEUT,
     $            IPNEUT,ABSPNEUT,PNEUT

C
C     Vector Parameter COMMON
C
      INTEGER*4   NPARVEC,IORGVEC(100),IPVEC(100),ICHVEC(100),
     $            IFLVEC(100)
      REAL*4      ABSPVEC(100),PVEC(3,100),POSV(3),AMASS(100)
      COMMON/COM2/NPARVEC,IPVEC,IORGVEC,ICHVEC,IFLVEC,
     $     ABSPVEC,PVEC,AMASS,POSV
C  
C     Flux info. common
C     

      integer*4 rlu
      External  rlu
      H(1)= rlu(i)

C
C --  INITIALIZATION
C
      CALL KZINIT
      CALL HLIMIT(-1*ISIZE)
      CALL SKOPENF(LUNI,0,'Z',IERR)

      IF (IARGC() .GE. 1) THEN
         CALL GETARG(1, FNAME_NTPL)
      ELSE
         FNAME_NTPL = 'genvec.nt'
      END IF

      CALL HROPEN(85,'CWNT',FNAME_NTPL,'N',LRECL,IERR)

      IF (IERR.ne.0) THEN
         WRITE(6,*) 'HROPEN:FAILED TO OPEN FILE ',FNAME_NTPL
         STOP
      ENDIF

      CALL HBNT(10,'NEUTVEC',' ')
      CALL HBNAME(10,'NEUTINFO',IEVENT,
     $     'NEV,NEUTMODE,NPNE[0,100],IPNE(NPNE),ABSPNE(NPNE),
     $PVECNE(3,NPNE)')

      CALL HBNAME(10,'VECINFO',NPARVEC,
     $     'NPVC[0,100],IPVC(NPVC),IORGVC(NPVC),ICRNVC(NPVC),
     $IFLGVC(NPVC),ABSPVC(NPVC),PVC(3,NPVC),AMASVC(NPVC),
     $POS(3)')

C      call hbname(10,'photons',nhits,'NHITS:I,POTOT:R')

C
C --  READ VECTOR
C

      DO 10 I=1,50000
C         CALL SKOPTN("31,30,29")
C         CALL SKREAD(LUNI,1000,2000,10,10)
         CALL KZREAD(LUNI,IERR)
         IF (IERR.eq.2) GOTO 1000
         IF (IERR.eq.1) GOTO 2000

         CALL KZBLOC('NEUT',LBANK)
         IF (LBANK.eq.0) GOTO 2000
         IF (IERR.EQ.2) GOTO 1000
         IF (IERR.EQ.1) GOTO 2000
         IF (MOD(I,1000).eq.0) then
             write(6,*)  ' *** ',i,' ***'
         endif
         IEVENT=I

         CALL NERDNEBK(POS)

         CALL VCCLCM
         CALL MCRDHD

         DO 60 N1=1,50
            IPNEUT(N1)=0
            ABSPNEUT(N1)=0.
            IPVEC(N1)=0
            IORGVEC(N1)=0
            ICHVEC(N1)=0
            ABSPVEC(N1)=0.
            AMASS(N1)=0.
            DO 70 N2=1,3
               PNEUT(N2,N1)=0.
               PVEC(N2,N1) =0.
 70         CONTINUE
 60      CONTINUE

         MODEN=MODENE         
         NPARNEUT = NUMNE
C         write(*,*) "MODEN=",MODEN," / NPAR=",NPAR

         DO 30 N1=1,NUMNE
            IPNEUT(N1)=IPNE(N1)
            ABSPNEUT(N1)=sqrt( PNE(1,N1)**2
     $                        +PNE(2,N1)**2
     $                        +PNE(3,N1)**2)
            
            PNEUT(1,N1)=PNE(1,N1)
            PNEUT(2,N1)=PNE(2,N1)
            PNEUT(3,N1)=PNE(3,N1)
            
C    write(*,*) "IP(",N,")=",IPNE(N),
C    $                 " / ABSP(",N,")=",ABSP(N),
C    $                 " / P(",N,")=(",PNE(1,N),",",PNE(2,N),",",
C    $                 PNE(3,N),")"
 30      CONTINUE
         CALL VCRDVCCM
         NPARVEC=NVC

         DO 40 N1=1,NVC
            IPVEC(N1)    = IPVC(N1)
            ICHVEC(N1)   = ICRNVC(N1)
            IORGVEC(N1)  = IORGVC(N1)
            IFLVEC(N1)   = IFLGVC(N1)
            
            ABSPVEC(N1)  = sqrt( PVC(1,N1)**2
     $                          +PVC(2,N1)**2
     $                          +PVC(3,N1)**2)
            
            PVEC(1,N1)   = PVC(1,N1)
            PVEC(2,N1)   = PVC(2,N1)
            PVEC(3,N1)   = PVC(3,N1)
            AMASS(N1)    = AMASVC(N1)

 40      CONTINUE

         posv(1)    = pvtxvc(1,1)
         posv(2)    = pvtxvc(2,1)
         posv(3)    = pvtxvc(3,1)

 100     continue

         CALL HFNT(10)
         CALL KZECLR
 10   CONTINUE
      GOTO 3000

 1000 WRITE(6,*) ' READ ERROR : ',LUNI
      GOTO 3000

 2000 WRITE(6,*) ' END OF FILE : ',LUNI

 3000 CONTINUE

      CALL HROUT(10,ICYCLE,' ')
      CALL HREND('CWNT')

C      CALL SKCLOSEF(LUNI)

      STOP
      END

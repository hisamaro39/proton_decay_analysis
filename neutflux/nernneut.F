************************************************************************
*     ------------------------------------
      SUBROUTINE NERNNEUT(ENEU, DIR, IPAR)
*     ------------------------------------
*
*     (Purpose)
*       Atmospehric neutrino event rate / kton yr
*
*     (Input)
*       ENEU   : NEUTRINO ENERGY ( GEV )
*       IPAR   : NEUTRINO TYPE
*                  12 : NEU-E
*                 -12 : NEU-E-BAR            
*                  14 : NEU-MU
*                 -14 : NEU-MU-BAR            
*
*     (Output)
*       DIR    : NEUTRINO DIRECTION
*
*     (Creation Date and Author)
*  --- OLD -----------------------------------------------------------
*  !   CONDITION :                                                   !
*  !   1. CALCULATION BY T.K.GAISSER.(JAN-84 VERSION)                !
*  !   2. SOLAR MAXIMUM.                                             !
*  !   (REF :  PRIVATE COMMUNICATION.)                               !
*  !                                                                 !
*  ! --MODIFICATION                                                  !
*  !   FOR ORIGINAL, SEE RNNEUT IN LCALP.                            !
*  !   AVERAGE OF SOLAR MAXIMUM AND MINIMUM IS INSTALLED BY M.TAKITA !
*  !       16-OCT-1987                                               !
*  --- OLD -----------------------------------------------------------
*       1988.11.24 ; T.KAJITA  LATEST GAISSER-FLUX IS INSTALLED (<3GEV)
*                              VOLKOVA FLUX IS INSTALLED FOR E(NEU)>10GEV
*                              SOLAR MIN,MAX OR ANY COMBINATION OF MIN/MAX
*                              IS AVAILABLE
*                                IF(SOLACT=1.) -- SOLAR MAX.
*                                IF(SOLACT=0.) -- SOLAR MIN.
*                               (IF(SOLACT=0.5) -- AVERAGE OF SOLAR MIN/MAX )
*                               NEU/(NEU-BAR) RATIO IS INSTALLED ACCORDING TO
*                               THE GAISSER'S TABLE
*                               ALSO NEU/(NEU-BAR) AT HIGH ENERGY IS ESTIMATED
*                               FROM THE SAME TABLE
*                               CORRECTION OF NEU(E)/NEU(MU) RATIO AT LOW E
*                               (<100MEV) IS ESTIMATED FROM THE HIGHER ENERGY
*                               CALCULATION
*       1994.02.01 ; K. KANEYUKI  FOR S.K.
*                    RNNEUT -> NERNNEUT
*       2007.08.27 ; G.Mitsuka Honda06 is supported     
*
************************************************************************
      DIMENSION DIR(3), ENELOW(10), PCOSL(11), ENEHI(3), PCOSH(8),
     &         FEMIN(11,10), FEMAX(11,10), FEBMIN(11,10), FEBMAX(11,10),
     &         FMMIN(11,10), FMMAX(11,10), FMBMIN(11,10), FMBMAX(11,10),
     &         FEHI(8,3),  FMHI(8,3)
#include "nework.h"
#include "nesolact.h"
#include "necardap.h"
      
C      COMMON /NESOLACT/ ISOL, SOLACT
C      COMMON/MODGEN/MODGE1,MODGE2,THRTOT,EDUM,DIRNU(3)
C
      DATA ENELOW /0.3, 0.5, 0.7, 0.9, 1.1, 1.3, 1.5, 1.7, 1.9, 3.0/
      DATA ENEHI  / 10., 30., 100. /
      DATA PCOSL / 1.0 ,  0.88,  0.63,  0.38, 0.13, 0.0 ,
     &            -0.13, -0.38, -0.63, -0.88, -1.0       /
      DATA PCOSH /1.0, 0.60, 0.40, 0.30, 0.20, 0.10, 0.05, 0.0 /
C
C NEU(E) FLUX AT SOLAR MIN. AND AT KAM. (T.K.GAISSER)
C                                               /M**2/SEC/SR/200MEV
      DATA FEMIN/.153E+03,.158E+03,.168E+03,.175E+03,.178E+03,
     &      .191E+03,.215E+03,.241E+03,.275E+03,.234E+03,.214E+03,
     &           .591E+02,.622E+02,.687E+02,.744E+02,.789E+02,
     &      .848E+02,.889E+02,.898E+02,.952E+02,.808E+02,.739E+02,
     &           .282E+02,.304E+02,.350E+02,.391E+02,.428E+02,
     &      .459E+02,.461E+02,.442E+02,.446E+02,.368E+02,.331E+02,
     &           .154E+02,.169E+02,.201E+02,.232E+02,.252E+02,
     &      .267E+02,.265E+02,.252E+02,.243E+02,.196E+02,.173E+02,
     &           .093E+02,.103E+02,.123E+02,.146E+02,.169E+02,
     &      .184E+02,.175E+02,.153E+02,.143E+02,.117E+02,.105E+02,
     &           .625E+01,.686E+01,.814E+01,.975E+01,.114E+02,
     &      .124E+02,.117E+02,.100E+02,.908E+01,.752E+01,.677E+01,
     &           .394E+01,.442E+01,.541E+01,.658E+01,.832E+01,
     &      .927E+01,.840E+01,.664E+01,.591E+01,.479E+01,.425E+01,
     &           .278E+01,.316E+01,.396E+01,.506E+01,.597E+01,
     &      .645E+01,.597E+01,.504E+01,.421E+01,.334E+01,.292E+01,
     &           .201E+01,.230E+01,.291E+01,.372E+01,.452E+01,
     &      .493E+01,.450E+01,.368E+01,.306E+01,.240E+01,.208E+01,
     &           .622E+00,.685E+00,.816E+00,1.07E+00,1.63E+00,
     &      1.89E+00,1.58E+00,1.03E+00,.811E+00,.686E+00,.626E+00 /
C
C NEU(E) FLUX AT SOLAR MAX. AND AT KAM. (T.K.GAISSER)
C                                               /M**2/SEC/SR/200MEV
      DATA FEMAX/.152E+03,.157E+03,.167E+03,.173E+03,.176E+03,
     &      .183E+03,.194E+03,.203E+03,.224E+03,.198E+03,.186E+03,
     &           .588E+02,.619E+02,.683E+02,.740E+02,.783E+02,
     &      .823E+02,.829E+02,.805E+02,.832E+02,.723E+02,.671E+02,
     &           .281E+02,.303E+02,.348E+02,.389E+02,.426E+02,
     &      .449E+02,.438E+02,.408E+02,.403E+02,.339E+02,.308E+02,
     &           .154E+02,.169E+02,.201E+02,.231E+02,.250E+02,
     &      .261E+02,.254E+02,.237E+02,.224E+02,.184E+02,.165E+02,
     &           .093E+02,.103E+02,.123E+02,.146E+02,.169E+02,
     &      .181E+02,.169E+02,.146E+02,.134E+02,.111E+02,.100E+02,
     &           .623E+01,.684E+01,.812E+01,.973E+01,.114E+02,
     &      .122E+02,.113E+02,.965E+01,.864E+01,.720E+01,.651E+01,
     &           .393E+01,.441E+01,.540E+01,.658E+01,.830E+01,
     &      .916E+01,.821E+01,.646E+01,.566E+01,.462E+01,.412E+01,
     &           .278E+01,.316E+01,.396E+01,.506E+01,.596E+01,
     &      .639E+01,.586E+01,.493E+01,.408E+01,.325E+01,.285E+01,
     &           .201E+01,.230E+01,.291E+01,.372E+01,.451E+01,
     &      .489E+01,.443E+01,.362E+01,.298E+01,.235E+01,.205E+01,
     &           .622E+00,.685E+00,.816E+00,1.07E+00,1.63E+00,
     &      1.89E+00,1.57E+00,1.02E+00,.798E+00,.676E+00,.617E+00  /
C
C
C NEU(E)BAR FLUX AT SOLAR MIN AND AT KAM. (T.K.GAISSER)
C                                               /M**2/SEC/SR/200MEV
      DATA FEBMIN/.146E+03,.150E+03,.159E+03,.166E+03,.170E+03,
     &       .181E+03,.198E+03,.214E+03,.239E+03,.207E+03,.192E+03,
     &            .545E+02,.573E+02,.631E+02,.687E+02,.732E+02,
     &       .781E+02,.799E+02,.786E+02,.821E+02,.706E+02,.651E+02,
     &            .253E+02,.272E+02,.312E+02,.356E+02,.385E+02,
     &       .408E+02,.405E+02,.384E+02,.376E+02,.316E+02,.287E+02,
     &            .139E+02,.151E+02,.177E+02,.203E+02,.229E+02,
     &       .244E+02,.235E+02,.213E+02,.201E+02,.167E+02,.151E+02,
     &            .812E+01,.889E+01,.105E+02,.125E+02,.145E+02,
     &       .156E+02,.146E+02,.127E+02,.115E+02,.960E+01,.869E+01,
     &            .501E+01,.562E+01,.689E+01,.833E+01,.958E+01,
     &       .102E+02,.960E+01,.836E+01,.736E+01,.599E+01,.533E+01,
     &            .327E+01,.371E+01,.462E+01,.558E+01,.700E+01,
     &       .772E+01,.698E+01,.558E+01,.488E+01,.392E+01,.346E+01,
     &            .245E+01,.275E+01,.338E+01,.415E+01,.497E+01,
     &       .539E+01,.495E+01,.412E+01,.353E+01,.286E+01,.254E+01,
     &            .181E+01,.203E+01,.248E+01,.307E+01,.386E+01,
     &       .424E+01,.380E+01,.301E+01,.256E+01,.210E+01,.188E+01,
     &            .446E+00,.516E+00,.662E+00,.877E+00,1.30E+00,
     &       1.51E+00,1.28E+00,.853E+00,.650E+00,.533E+00,.477E+00  /
C
C NEU(E)BAR FLUX AT SOLAR MAX AND AT KAM. (T.K.GAISSER)
C                                               /M**2/SEC/SR/200MEV
      DATA FEBMAX/.145E+03,.149E+03,.158E+03,.165E+03,.168E+03,
     &       .175E+03,.182E+03,.186E+03,.202E+03,.180E+03,.169E+03,
     &            .542E+02,.570E+02,.628E+02,.683E+02,.727E+02,
     &       .762E+02,.755E+02,.720E+02,.735E+02,.645E+02,.602E+02,
     &            .252E+02,.271E+02,.311E+02,.355E+02,.382E+02,
     &       .400E+02,.389E+02,.361E+02,.347E+02,.295E+02,.270E+02,
     &            .138E+02,.150E+02,.176E+02,.202E+02,.228E+02,
     &       .241E+02,.228E+02,.203E+02,.190E+02,.159E+02,.144E+02,
     &            .809E+01,.887E+01,.105E+02,.124E+02,.144E+02,
     &       .154E+02,.143E+02,.123E+02,.110E+02,.926E+01,.842E+01,
     &            .500E+01,.561E+01,.688E+01,.832E+01,.955E+01,
     &       .101E+02,.943E+01,.816E+01,.712E+01,.581E+01,.518E+01,
     &            .327E+01,.371E+01,.462E+01,.557E+01,.699E+01,
     &       .767E+01,.688E+01,.547E+01,.474E+01,.382E+01,.338E+01,
     &            .245E+01,.275E+01,.338E+01,.414E+01,.496E+01,
     &       .535E+01,.488E+01,.405E+01,.345E+01,.280E+01,.249E+01,
     &            .182E+01,.203E+01,.247E+01,.307E+01,.385E+01,
     &       .422E+01,.377E+01,.298E+01,.251E+01,.206E+01,.184E+01,
     &            .446E+00,.516E+00,.662E+00,.877E+00,1.30E+00,
     &       1.50E+00,1.27E+00,.848E+00,.644E+00,.528E+00,.472E+00  /
C
C
C NEU(MU) FLUX AT SOLAR MIN. AND AT KAM. (T.K.GAISSER)
C                                               /M**2/SEC/SR/200MEV
      DATA FMMIN/.334E+03,.338E+03,.347E+03,.349E+03,.353E+03,
     &      .374E+03,.415E+03,.459E+03,.535E+03,.479E+03,.452E+03,
     &           .135E+03,.137E+03,.142E+03,.147E+03,.150E+03,
     &      .157E+03,.165E+03,.169E+03,.185E+03,.169E+03,.161E+03,
     &           .673E+02,.690E+02,.726E+02,.752E+02,.779E+02,
     &      .809E+02,.821E+02,.812E+02,.872E+02,.800E+02,.765E+02,
     &           .380E+02,.390E+02,.411E+02,.442E+02,.459E+02,
     &      .475E+02,.474E+02,.460E+02,.469E+02,.433E+02,.416E+02,
     &           .233E+02,.240E+02,.255E+02,.278E+02,.295E+02,
     &      .305E+02,.299E+02,.284E+02,.281E+02,.260E+02,.250E+02,
     &           .158E+02,.163E+02,.174E+02,.186E+02,.201E+02,
     &      .209E+02,.201E+02,.187E+02,.186E+02,.172E+02,.165E+02,
     &           .108E+02,.112E+02,.120E+02,.130E+02,.145E+02,
     &      .153E+02,.144E+02,.128E+02,.125E+02,.116E+02,.112E+02,
     &           .752E+01,.787E+01,.859E+01,.929E+01,.106E+02,
     &      .112E+02,.105E+02,.918E+01,.883E+01,.805E+01,.768E+01,
     &           .545E+01,.572E+01,.629E+01,.717E+01,.785E+01,
     &      .814E+01,.771E+01,.701E+01,.638E+01,.579E+01,.551E+01,
     &           1.89E+00,1.95E+00,2.07E+00,2.25E+00,2.70E+00,
     &      2.89E+00,2.62E+00,2.20E+00,2.04E+00,1.96E+00,1.92E+00  /
C
C NEU(MU) FLUX AT SOLAR MAX. AND AT KAM. (T.K.GAISSER)
C                                               /M**2/SEC/SR/200MEV
      DATA FMMAX/.332E+03,.336E+03,.344E+03,.347E+03,.349E+03,
     &      .360E+03,.379E+03,.396E+03,.448E+03,.413E+03,.396E+03,
     &           .134E+03,.136E+03,.141E+03,.146E+03,.149E+03,
     &      .154E+03,.156E+03,.154E+03,.166E+03,.154E+03,.148E+03,
     &           .671E+02,.688E+02,.724E+02,.750E+02,.774E+02,
     &      .792E+02,.787E+02,.765E+02,.807E+02,.750E+02,.723E+02,
     &           .379E+02,.389E+02,.410E+02,.441E+02,.457E+02,
     &      .467E+02,.459E+02,.440E+02,.442E+02,.413E+02,.399E+02,
     &           .231E+02,.239E+02,.255E+02,.278E+02,.294E+02,
     &      .302E+02,.292E+02,.275E+02,.269E+02,.250E+02,.241E+02,
     &           .158E+02,.163E+02,.174E+02,.186E+02,.200E+02,
     &      .206E+02,.197E+02,.182E+02,.180E+02,.168E+02,.162E+02,
     &           .108E+02,.112E+02,.120E+02,.130E+02,.145E+02,
     &      .152E+02,.142E+02,.126E+02,.122E+02,.114E+02,.110E+02,
     &           .753E+01,.787E+01,.858E+01,.928E+01,.106E+02,
     &      .112E+02,.104E+02,.905E+01,.867E+01,.794E+01,.759E+01,
     &           .545E+01,.572E+01,.629E+01,.716E+01,.785E+01,
     &      .812E+01,.766E+01,.694E+01,.630E+01,.574E+01,.547E+01,
     &           1.89E+00,1.95E+00,2.07E+00,2.25E+00,2.70E+00,
     &      2.89E+00,2.62E+00,2.19E+00,2.03E+00,1.95E+00,1.91E+00  /
C
C NEU(MU)BAR FLUX AT SOLAR MIN AND AT KAM (T.K.GAISSER)
C                                               /M**2/SEC/SR/200MEV
      DATA FMBMIN/.333E+03,.337E+03,.345E+03,.350E+03,.352E+03,
     &       .371E+03,.414E+03,.461E+03,.534E+03,.477E+03,.450E+03,
     &            .132E+03,.135E+03,.141E+03,.146E+03,.149E+03,
     &       .155E+03,.163E+03,.168E+03,.185E+03,.167E+03,.158E+03,
     &            .654E+02,.672E+02,.710E+02,.748E+02,.775E+02,
     &       .809E+02,.823E+02,.812E+02,.856E+02,.777E+02,.739E+02,
     &            .366E+02,.380E+02,.409E+02,.434E+02,.455E+02,
     &       .473E+02,.472E+02,.456E+02,.467E+02,.423E+02,.402E+02,
     &            .226E+02,.236E+02,.256E+02,.274E+02,.291E+02,
     &       .302E+02,.297E+02,.282E+02,.282E+02,.254E+02,.241E+02,
     &            .147E+02,.154E+02,.168E+02,.185E+02,.202E+02,
     &       .211E+02,.202E+02,.186E+02,.180E+02,.164E+02,.156E+02,
     &            .102E+02,.106E+02,.115E+02,.126E+02,.144E+02,
     &       .153E+02,.144E+02,.126E+02,.122E+02,.110E+02,.104E+02,
     &            .702E+01,.745E+01,.834E+01,.918E+01,.104E+02,
     &       .110E+02,.103E+02,.909E+01,.858E+01,.764E+01,.719E+01,
     &            .516E+01,.547E+01,.612E+01,.694E+01,.800E+01,
     &       .849E+01,.788E+01,.683E+01,.628E+01,.559E+01,.526E+01,
     &            1.67E+00,1.77E+00,1.97E+00,2.22E+00,2.76E+00,
     &       3.01E+00,2.69E+00,2.15E+00,1.95E+00,1.76E+00,1.67E+00  /
C
C NEU(MU)BAR FLUX AT SOLAR MAX AND AT KAM (T.K.GAISSER)
C                                               /M**2/SEC/SR/200MEV
      DATA FMBMAX/.330E+03,.334E+03,.343E+03,.347E+03,.348E+03,
     &       .359E+03,.379E+03,.396E+03,.446E+03,.411E+03,.394E+03,
     &            .133E+03,.135E+03,.140E+03,.145E+03,.147E+03,
     &       .151E+03,.154E+03,.154E+03,.165E+03,.153E+03,.147E+03,
     &            .652E+02,.670E+02,.707E+02,.745E+02,.771E+02,
     &       .792E+02,.787E+02,.763E+02,.790E+02,.729E+02,.700E+02,
     &            .365E+02,.379E+02,.408E+02,.432E+02,.453E+02,
     &       .465E+02,.456E+02,.435E+02,.440E+02,.403E+02,.385E+02,
     &            .225E+02,.235E+02,.255E+02,.274E+02,.290E+02,
     &       .298E+02,.289E+02,.272E+02,.269E+02,.245E+02,.233E+02,
     &            .147E+02,.154E+02,.168E+02,.184E+02,.201E+02,
     &       .208E+02,.198E+02,.181E+02,.174E+02,.159E+02,.152E+02,
     &            .100E+02,.105E+02,.115E+02,.126E+02,.143E+02,
     &       .151E+02,.141E+02,.123E+02,.119E+02,.108E+02,.103E+02,
     &            .702E+01,.745E+01,.834E+01,.918E+01,.104E+02,
     &       .109E+02,.102E+02,.895E+01,.842E+01,.753E+01,.710E+01,
     &            .516E+01,.547E+01,.611E+01,.693E+01,.799E+01,
     &       .844E+01,.780E+01,.675E+01,.617E+01,.551E+01,.519E+01,
     &            1.67E+00,1.77E+00,1.97E+00,2.22E+00,2.76E+00,
     &       3.00E+00,2.68E+00,2.14E+00,1.94E+00,1.75E+00,1.66E+00  /
C
C NEU(E)+NEU(E)-BAR FLUX AT HIGH ENERGIES (VOLKOVA) /M**2/SEC/SR/GEV
      DATA FEHI/ .681E-1,  1.17E-1,  1.69E-1,  2.13E-1,  2.83E-1,
     &           3.66E-1,  3.99E-1,  4.18E-1,
     &           1.45E-3,  2.49E-3,  3.74E-3,  4.90E-3,  7.00E-3,
     &           10.1E-3,  11.6E-3,  12.7E-3,
     &           2.15E-5,  3.40E-5,  5.00E-5,  6.70E-5,  10.1E-5,
     &           16.5E-5,  20.6E-5,  24.2E-5          /
C
C NEU(MU)+NEU(MU)-BAR FLUX AT HIGH ENERGIES (VOLKOVA) /M**2/SEC/SR/GEV
      DATA FMHI/ 4.38E-1,  5.19E-1,  5.91E-1,  6.50E-1,  7.33E-1,
     &           8.28E-1,  8.65E-1,  8.86E-1,
     &           1.53E-2,  1.86E-2,  2.17E-2,  2.42E-2,  2.79E-2,
     &           3.25E-2,  3.44E-2,  3.57E-2,
     &           3.91E-4,  4.95E-4,  5.89E-4,  6.62E-4,  7.70E-4,
     &           9.19E-4,  9.91E-4,  10.5E-4           /
C
C --  STATEMENT FUNCTION
C             NEU/(NEU+NEUBAR) AT E>3GEV FOR NEU(E) AND NEU(MU)

      FNENEB(COSTH)=0.5575
      FNMNMB(COSTH)=0.5025+ 0.045*ABS(COSTH)

      IF(IPAR.EQ.16) THEN
         IPAR_FLUX=14
      ELSEIF(IPAR.EQ.-16) THEN
         IPAR_FLUX=-14
      ELSE
         IPAR_FLUX=IPAR
      ENDIF
C
C-- CHECK KIND OF FLUX
C
      IF (NEFLXFLG.eq.0.or.NEFLXFLG.eq.6) THEN
         CALL NERNNUDH(ENEU,DIR,IPAR_FLUX)
         RETURN
      ENDIF

      IF (NEFLXFLG.eq.1) THEN
         CALL NERNNUDG(ENEU,DIR,IPAR_FLUX)
         RETURN
      ENDIF

      IF (NEFLXFLG.eq.5) THEN
         CALL NERNNUDH_1D(ENEU,DIR,IPAR_FLUX)
         RETURN
      ENDIF

      IF (NEFLXFLG.eq.2.or.NEFLXFLG.eq.3.or.
     &     NEFLXFLG.eq.4) THEN
         WRITE(*,*) "NEFLXFLG =",NEFLXFLG," is not supported."
         STOP
      ENDIF
      
CC Y.H COMMENT OUT MODGEN RELATED
C
C     
C --  FIXED PARTICLE DIRECTION
C
C      IF ( MODGEN.EQ.-1 ) THEN
C         DIR(1)=DIRINP(1)
C         DIR(2)=DIRINP(2)
C         DIR(3)=DIRINP(3)
C         RETURN
C      ENDIF

C     
C --  CONSIDER SOLAR ACTIVITY
C      
      CALL VZERO(DIR,3)
      IF(ISOL.EQ.0) SOLACT=0.5
      IF(SOLACT.LT.0.) SOLACT=0.
      IF(SOLACT.GT.1.) SOLACT=1.

C
C --  DETERMINE ENERGY BINS
C
      DO 1 I=1,10
         IMAX=I
         DENEU = ENEU-ENELOW(I)
         IF(DENEU.LT.0.) THEN
            ILOW=I-1
            IF(ILOW.EQ.0) ILOW=1
            IHIGH=I
            IENERG=1
            GO TO 10
         END IF
 1    CONTINUE

      DO 2 J=1,3
         JMAX=J
         DENEU = ENEU-ENEHI(J)
         IF(DENEU.LT.0.) THEN
            ILOW=J-1
            IHIGH=J
            IENERG=3
            IF(ILOW.EQ.0) THEN
               ILOW=IMAX
               IENERG=2
            END IF
            GO TO 10
         END IF
 2    CONTINUE
      ILOW=JMAX
      IHIGH=JMAX
      IENERG=3
 10   CONTINUE
C
C
C      GO TO (100,200,300 ),IENERG
C
C
C -- E(NEU) LESS THAN 3 GEV
C
      IF (IENERG.EQ.1) THEN
 100     CONTINUE
         COST = 2.*RLU(DUM)-1.
         DO 110 ICOS=2,11
            IF (COST.GT.PCOSL(ICOS)) THEN
               ICOS1 = ICOS-1
               ICOS2 = ICOS
               DCOS1 = - ( COST-PCOSL(ICOS-1) )
               DCOS2 =     COST-PCOSL(ICOS)
               GO TO 111
            END IF
 110     CONTINUE
         RETURN
 111     CONTINUE

         IF (IPAR_FLUX .EQ. 12) THEN
            FLLOW = FEMIN(6,ILOW )
            FLHIGH= FEMIN(6,IHIGH)
            FLCOSL= SOLACT*( DCOS2*FEMAX(ICOS1,ILOW )
     &                       +DCOS1*FEMAX(ICOS2,ILOW ))/(DCOS1+DCOS2) +
     &          (1.-SOLACT)*( DCOS2*FEMIN(ICOS1,ILOW )
     &                       +DCOS1*FEMIN(ICOS2,ILOW ))/(DCOS1+DCOS2)
            FLCOSH= SOLACT*( DCOS2*FEMAX(ICOS1,IHIGH)
     &                       +DCOS1*FEMAX(ICOS2,IHIGH))/(DCOS1+DCOS2) +
     &          (1.-SOLACT)*( DCOS2*FEMIN(ICOS1,IHIGH)
     &                       +DCOS1*FEMIN(ICOS2,IHIGH))/(DCOS1+DCOS2)
         END IF
         IF(IPAR_FLUX.EQ.-12) THEN
            FLLOW = FEBMIN(6,ILOW )
            FLHIGH= FEBMIN(6,IHIGH)
            FLCOSL= SOLACT*( DCOS2*FEBMAX(ICOS1,ILOW )
     &                       +DCOS1*FEBMAX(ICOS2,ILOW ))/(DCOS1+DCOS2) +
     &          (1.-SOLACT)*( DCOS2*FEBMIN(ICOS1,ILOW )
     &                       +DCOS1*FEBMIN(ICOS2,ILOW ))/(DCOS1+DCOS2)
            FLCOSH= SOLACT*( DCOS2*FEBMAX(ICOS1,IHIGH)
     &                       +DCOS1*FEBMAX(ICOS2,IHIGH))/(DCOS1+DCOS2) +
     &          (1.-SOLACT)*( DCOS2*FEBMIN(ICOS1,IHIGH)
     &                       +DCOS1*FEBMIN(ICOS2,IHIGH))/(DCOS1+DCOS2)
         END IF
         IF(IPAR_FLUX.EQ. 14) THEN
            FLLOW = FMMIN(6,ILOW )
            FLHIGH= FMMIN(6,IHIGH)
            FLCOSL= SOLACT*( DCOS2*FMMAX(ICOS1,ILOW )
     &                       +DCOS1*FMMAX(ICOS2,ILOW ))/(DCOS1+DCOS2) +
     &          (1.-SOLACT)*( DCOS2*FMMIN(ICOS1,ILOW )
     &                       +DCOS1*FMMIN(ICOS2,ILOW ))/(DCOS1+DCOS2)
            FLCOSH= SOLACT*( DCOS2*FMMAX(ICOS1,IHIGH)
     &                       +DCOS1*FMMAX(ICOS2,IHIGH))/(DCOS1+DCOS2) +
     &          (1.-SOLACT)*( DCOS2*FMMIN(ICOS1,IHIGH)
     &                       +DCOS1*FMMIN(ICOS2,IHIGH))/(DCOS1+DCOS2)
         END IF
         IF(IPAR_FLUX.EQ.-14) THEN
            FLLOW = FMBMIN(6,ILOW )
            FLHIGH= FMBMIN(6,IHIGH)
            FLCOSL= SOLACT*( DCOS2*FMBMAX(ICOS1,ILOW )
     &                       +DCOS1*FMBMAX(ICOS2,ILOW ))/(DCOS1+DCOS2) +
     &          (1.-SOLACT)*( DCOS2*FMBMIN(ICOS1,ILOW )
     &                       +DCOS1*FMBMIN(ICOS2,ILOW ))/(DCOS1+DCOS2)
            FLCOSH= SOLACT*( DCOS2*FMBMAX(ICOS1,IHIGH)
     &                       +DCOS1*FMBMAX(ICOS2,IHIGH))/(DCOS1+DCOS2) +
     &          (1.-SOLACT)*( DCOS2*FMBMIN(ICOS1,IHIGH)
     &                       +DCOS1*FMBMIN(ICOS2,IHIGH))/(DCOS1+DCOS2)
         END IF
         IF(ILOW.NE.IHIGH) THEN
            BNORM = ALOG(FLLOW/FLHIGH)/ALOG(ENELOW(ILOW)/ENELOW(IHIGH))
            ANORM = FLLOW / ENELOW(ILOW)**BNORM
            FXNORM = ANORM*ENEU**BNORM *2.
            B = ALOG(FLCOSL/FLCOSH)/ALOG(ENELOW(ILOW)/ENELOW(IHIGH))
            A = FLCOSL / ENELOW(ILOW)**B
            FLUXNE= A*ENEU**B
         ELSE
            FXNORM=FLLOW*2.
            FLUXNE=FLCOSL
         END IF
         XRN=RLU(DUM)
         IF(XRN.GT.(FLUXNE/FXNORM)) GO TO 100

C
C -- E(NEU) BETWEEN 3 GEV AND 10 GEV
C
      ELSE IF (IENERG.EQ.2) THEN 
 200     CONTINUE
         COST = 2.*RLU(DUM)-1.
         DO 210 ICOS=2,11
            IF (COST.GT.PCOSL(ICOS)) THEN
               ICOS1 = ICOS-1
               ICOS2 = ICOS
               DCOS1 = - ( COST-PCOSL(ICOS-1) )
               DCOS2 =     COST-PCOSL(ICOS)
               GO TO 211
            END IF
 210     CONTINUE
         RETURN
 211     CONTINUE

         DO 220 ICOS=2,8
            IF (ABS(COST).GT.PCOSH(ICOS)) THEN
               ICOS3 = ICOS-1
               ICOS4 = ICOS
               DCOS3 = - ( ABS(COST)-PCOSH(ICOS-1) )
               DCOS4 =     ABS(COST)-PCOSH(ICOS)
               GO TO 221
            END IF
 220     CONTINUE
         RETURN
 221     CONTINUE

         IF(IPAR_FLUX.EQ. 12) THEN
            FLLOW = FEMIN(6,ILOW )
            FLHIGH= FEHI(8,IHIGH)*FNENEB(0.)/5.
            FLCOSL= SOLACT*( DCOS2*FEMAX(ICOS1,ILOW )
     &                       +DCOS1*FEMAX(ICOS2,ILOW ))/(DCOS1+DCOS2) +
     &          (1.-SOLACT)*( DCOS2*FEMIN(ICOS1,ILOW )
     &                       +DCOS1*FEMIN(ICOS2,ILOW ))/(DCOS1+DCOS2)
            FLCOSH= ( DCOS4*FEHI(ICOS3,IHIGH)+DCOS3*FEHI(ICOS4,IHIGH))
     &                      /(DCOS3+DCOS4)*FNENEB(COST)/5.
         END IF
         IF(IPAR_FLUX.EQ.-12) THEN
            FLLOW = FEBMIN(6,ILOW )
            FLHIGH= FEHI(8,IHIGH)*(1.-FNENEB(0.))/5.
            FLCOSL= SOLACT*( DCOS2*FEBMAX(ICOS1,ILOW )
     &                       +DCOS1*FEBMAX(ICOS2,ILOW ))/(DCOS1+DCOS2) +
     &          (1.-SOLACT)*( DCOS2*FEBMIN(ICOS1,ILOW )
     &                       +DCOS1*FEBMIN(ICOS2,ILOW ))/(DCOS1+DCOS2)
            FLCOSH= ( DCOS4*FEHI(ICOS3,IHIGH)+DCOS3*FEHI(ICOS4,IHIGH))
     &                      /(DCOS3+DCOS4)*(1.-FNENEB(COST))/5.
         END IF
         IF(IPAR_FLUX.EQ. 14) THEN
            FLLOW = FMMIN(6,ILOW )
            FLHIGH= FMHI(8,IHIGH)*FNMNMB(0.)/5.
            FLCOSL= SOLACT*( DCOS2*FMMAX(ICOS1,ILOW )
     &                       +DCOS1*FMMAX(ICOS2,ILOW ))/(DCOS1+DCOS2) +
     &          (1.-SOLACT)*( DCOS2*FMMIN(ICOS1,ILOW )
     &                       +DCOS1*FMMIN(ICOS2,ILOW ))/(DCOS1+DCOS2)
            FLCOSH= ( DCOS4*FMHI(ICOS3,IHIGH)+DCOS3*FMHI(ICOS4,IHIGH))
     &                      /(DCOS3+DCOS4)*FNMNMB(COST)/5.
         END IF
         IF(IPAR_FLUX.EQ.-14) THEN
            FLLOW = FMBMIN(6,ILOW )
            FLHIGH= FMHI(8,IHIGH)*(1.-FNMNMB(0.))/5.
            FLCOSL= SOLACT*( DCOS2*FMBMAX(ICOS1,ILOW )
     &                       +DCOS1*FMBMAX(ICOS2,ILOW ))/(DCOS1+DCOS2) +
     &          (1.-SOLACT)*( DCOS2*FMBMIN(ICOS1,ILOW )
     &                       +DCOS1*FMBMIN(ICOS2,ILOW ))/(DCOS1+DCOS2)
            FLCOSH= ( DCOS4*FMHI(ICOS3,IHIGH)+DCOS3*FMHI(ICOS4,IHIGH))
     &                      /(DCOS3+DCOS4)*(1.-FNMNMB(COST))/5.
         END IF
         BNORM = ALOG(FLLOW/FLHIGH)/ALOG(ENELOW(ILOW)/ENEHI(IHIGH))
         ANORM = FLLOW / ENELOW(ILOW)**BNORM
         FXNORM = ANORM*ENEU**BNORM *2.
         B = ALOG (FLCOSL/FLCOSH) / ALOG (ENELOW(ILOW)/ENEHI(IHIGH))
         A = FLCOSL / ENELOW(ILOW)**B
         FLUXNE= A*ENEU**B
         XRN=RLU(DUM)
         IF(XRN.GT.(FLUXNE/FXNORM)) GO TO 200
C
C -- E(NEU) LARGER THAN 10 GEV
C
      ELSE IF (IENERG.EQ.3) THEN
 300     CONTINUE
         COST = 2.*RLU(DUM)-1.
         DO 320 ICOS=2,8
            IF (ABS(COST).GT.PCOSH(ICOS)) THEN
               ICOS3 = ICOS-1
               ICOS4 = ICOS
               DCOS3 = - ( ABS(COST)-PCOSH(ICOS-1) )
               DCOS4 =     ABS(COST)-PCOSH(ICOS)
               GO TO 321
            END IF
 320     CONTINUE
         RETURN
 321     CONTINUE

         IF(IPAR_FLUX.EQ. 12) THEN
            FLLOW = FEHI(8,ILOW )*FNENEB(0.)/5.
            FLHIGH= FEHI(8,IHIGH)*FNENEB(0.)/5.
            FLCOSL= ( DCOS4*FEHI(ICOS3,ILOW )+DCOS3*FEHI(ICOS4,ILOW ))
     &                      /(DCOS3+DCOS4)*FNENEB(COST)/5.
            FLCOSH= ( DCOS4*FEHI(ICOS3,IHIGH)+DCOS3*FEHI(ICOS4,IHIGH))
     &                      /(DCOS3+DCOS4)*FNENEB(COST)/5.
         END IF
         IF(IPAR_FLUX.EQ.-12) THEN
            FLLOW = FEHI(8,ILOW )*(1.-FNENEB(0.))/5.
            FLHIGH= FEHI(8,IHIGH)*(1.-FNENEB(0.))/5.
            FLCOSL= ( DCOS4*FEHI(ICOS3,ILOW )+DCOS3*FEHI(ICOS4,ILOW ))
     &                      /(DCOS3+DCOS4)*(1.-FNENEB(COST))/5.
            FLCOSH= ( DCOS4*FEHI(ICOS3,IHIGH)+DCOS3*FEHI(ICOS4,IHIGH))
     &                      /(DCOS3+DCOS4)*(1.-FNENEB(COST))/5.
         END IF
         IF(IPAR_FLUX.EQ. 14) THEN
            FLLOW = FMHI(8,ILOW )*FNMNMB(0.)/5.
            FLHIGH= FMHI(8,IHIGH)*FNMNMB(0.)/5.
            FLCOSL= ( DCOS4*FMHI(ICOS3,ILOW )+DCOS3*FMHI(ICOS4,ILOW ))
     &                      /(DCOS3+DCOS4)*FNMNMB(COST)/5.
            FLCOSH= ( DCOS4*FMHI(ICOS3,IHIGH)+DCOS3*FMHI(ICOS4,IHIGH))
     &                      /(DCOS3+DCOS4)*FNMNMB(COST)/5.
         END IF
         IF(IPAR_FLUX.EQ.-14) THEN
            FLLOW = FMHI(8,ILOW )*(1.-FNMNMB(0.))/5.
            FLHIGH= FMHI(8,IHIGH)*(1.-FNMNMB(0.))/5.
            FLCOSL= ( DCOS4*FMHI(ICOS3,ILOW )+DCOS3*FMHI(ICOS4,ILOW ))
     &                      /(DCOS3+DCOS4)*(1.-FNMNMB(COST))/5.
            FLCOSH= ( DCOS4*FMHI(ICOS3,IHIGH)+DCOS3*FMHI(ICOS4,IHIGH))
     &                      /(DCOS3+DCOS4)*(1.-FNMNMB(COST))/5.
         END IF
         IF(ILOW.NE.IHIGH) THEN
            BNORM = ALOG(FLLOW/FLHIGH)/ALOG(ENEHI(ILOW)/ENEHI(IHIGH))
            ANORM = FLLOW / ENEHI(ILOW)**BNORM
            FXNORM = ANORM*ENEU**BNORM *2.
            B = ALOG (FLCOSL/FLCOSH) / ALOG (ENEHI(ILOW)/ENEHI(IHIGH))
            A = FLCOSL / ENEHI(ILOW)**B
            FLUXNE= A*ENEU**B
         ELSE
            FXNORM=FLLOW*2.
            FLUXNE=FLCOSL
         END IF
         XRN=RLU(DUM)
         IF(XRN.GT.(FLUXNE/FXNORM)) GO TO 300
      ENDIF
C
C -- NOW SET NEUTRINO DIRECTION
C
      SINT=SQRT(1.-COST**2)
      PHI=RLU(DUM)*6.283184
      DIR(1)=COS(PHI)*SINT
      DIR(2)=SIN(PHI)*SINT
      DIR(3)=-COST
      RETURN
      END

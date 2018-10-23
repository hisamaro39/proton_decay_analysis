************************************************************************
*     --------
*     NECARDAP.H
*     --------
*
*     (Purpose)
*       COMMON for CARD on NEUT for atmospheric neutrino
*
*     (Variables)
**
*       ITAUFLG    : control Tau run mode
*       NEWATERVOL : Set Simulated water volume
*                    0(Full volume), 1(ID-Dwall), 2(OD)
*     
*     (Creation Date and Author)
*	2007.04.24 ; G.Mitsuka necard -> necardap(for ATMPD)
*		     add selection of coherent pion model
*       2007.08.?? ; G.Mitsuka deleted overlaps with necard.F
*       2007.12.05 ; G.Mitsuka added NEWATERVOL
************************************************************************
      INTEGER NEFLXFLG,ITAUFLG,NEWATERVOL
      REAL    NEYEARS

      COMMON/NEUTCARDAP/NEFLXFLG,NEYEARS,ITAUFLG,NEWATERVOL

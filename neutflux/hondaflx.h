******************************************
*        HFLXL: READ table LE
*        HFLXM: READ table HE
*        HFLX:
*        HFLXA: 
******************************************
      REAL*4 HFLXL(72,20,4,2),HFLXM(107,20,12,4,2),
     $       HFLX(141,4,2),HFLXA(179,4,2)
      REAL*4 HONE(141),HONEL(78),HONEM(107)
      REAL*4 HONEA(78+107)
      REAL*4 HFLXLNRM(4,2)

      COMMON /HONFX/HONE ,HFLX,
     $              HONEA,HFLXA,
     $              HONEL,HFLXL,
     $              HONEM,HFLXM,
     $              HFLXLNRM

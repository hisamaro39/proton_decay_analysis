      REAL*4 HFLXL_1D(72,20,4,2),HFLXM_1D(50,20,12,4,2),HFLXH_1D(31,20,4),
     $       HFLX_1D(100,4,2),HFLXA_1D(153,4,2)
      REAL*4 HONE_1D(100),HONEL_1D(72),HONEM_1D(50),HONEH_1D(31)
      REAL*4 HONEA_1D(72+50+31)
      REAL*4 HFLXLNRM_1D(4,2)
      COMMON /HONFX_1D/HONE_1D ,HFLX_1D,
     $              HONEA_1D,HFLXA_1D,
     $              HONEL_1D,HFLXL_1D,
     $              HONEM_1D,HFLXM_1D,
     $              HONEH_1D,HFLXH_1D,
     $              HFLXLNRM_1D





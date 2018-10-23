************************************************************
*
*      upmu.h
*
*     upmu_eth: upmu energy threshold(GeV)
*     dens_rck: rock density(g/cm3)
*     rrock:    radius of interaction point(cm)
*     upmu_diam:   diameter of tutu(cm)
*     upmu_rmin:     radius of holly region(cm)
*     upmu_maxene:  MAX energy of neutrino(GeV)
************************************************************
      real*4 upmu_eth,rrock,upmu_diam,upmu_rmin,dens_rck,pi
      real*4 upmu_maxene
      parameter( upmu_eth=1.6 )
      parameter( dens_rck=2.7 )
      parameter( rrock=400000. )
      parameter( upmu_diam=10000. )
      parameter( upmu_rmin=3000. )
      parameter( upmu_maxene=90000. )
      PARAMETER( PI=3.141592653 )

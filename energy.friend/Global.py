#!/usr/bin/python

import math


class Callable:
  def __init__( self, anycallable ):
    self.__call__ = anycallable


class OscTypes :

#MRO mode
     nTypes = 19 
#    nTypes = 18 

     SubGeV_elike_0dcy            =  1
     SubGeV_elike_1dcy            =  2
     SubGeV_SingleRing_pi0like    =  3
     SubGeV_mulike_0dcy           =  4
     SubGeV_mulike_1dcy           =  5
     SubGeV_mulike_2dcy           =  6
     SubGeV_pi0like               =  7
     MultiGeV_elike_nue           =  8
     MultiGeV_elike_nuebar        =  9
     MultiGeV_mulike              = 10
     MultiRing_elike_nue          = 11
     MultiRing_elike_nuebar       = 12
     MultiRing_mulike             = 13
     MultiRing_Other              = 14
     PCStop                       = 15
     PCThru                       = 16
     UpStop_mu                    = 17
     UpThruNonShower_mu           = 18 
     UpThruShower_mu              = 19 


     def GetZBin( cz , itype ):

        ncz      = 10.0 
        czstart  = -1.0
        czstop   =  1.0 
     
        if( cz >= 1.0 ) : cz = 0.9999999

        cmp  =  ( itype -1 ) % OscTypes.nTypes 
        cmp += 1 
 
        # print itype, OscTypes.nTypes,  cmp , OscTypes.SubGeV_SingleRing_pi0like 
        if ( cmp == OscTypes.SubGeV_SingleRing_pi0like  ) : ncz = 1.0 
        if ( cmp == OscTypes.SubGeV_elike_1dcy          ) : ncz = 1.0 
        if ( cmp == OscTypes.SubGeV_mulike_2dcy         ) : ncz = 1.0 
        if ( cmp == OscTypes.SubGeV_pi0like             ) : ncz = 1.0 

        if ( cmp == OscTypes.UpStop_mu           or \
             cmp == OscTypes.UpThruNonShower_mu  or \
             cmp == OscTypes.UpThruShower_mu     ) :
           czstop = 0.0  

        binWidth = (czstop - czstart)/ncz      
        bin      = int( (cz - czstart) /binWidth ) 

        return bin

     GetZBin = Callable( GetZBin )

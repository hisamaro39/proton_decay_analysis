#!/usr/bin/python

import math

class EAveContainer:

   def __init__( self , type, zbin, amom, pnu, ipnu , event ):
     self.type   = type
     self.amom   = amom
     self.pnu    = pnu
     self.ipnu   = ipnu 
     self.event  = event
     self.zbin   = zbin 

   def __repr__(self):
     return repr((self.type, self.zbin, self.amom, self.pnu, self.ipnu , self.event ))




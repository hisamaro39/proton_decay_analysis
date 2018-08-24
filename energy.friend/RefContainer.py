#!/usr/bin/python

import math

class RefContainer:

   def __init__( self , type, zbin, ipnu, dict ):
     self.type   = type
     self.zbin   = zbin 
     self.ipnu   = ipnu 
     self.dict   = dict 

   def __repr__(self):
     return repr((self.type, self.zbin, self.ipnu, self.dict ))


   

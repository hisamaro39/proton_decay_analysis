#!/usr/bin/python


from operator import itemgetter, attrgetter
import struct, array
import math
import sys
sys.argv.append("-b")
from ROOT import *
sys.argv.pop(-1)

from EAveContainer import EAveContainer
from RefContainer  import RefContainer
from Global        import OscTypes   


#
#  Routine for creating a friend tree of neighbor
#  neutrino energies based on type and reconstructed
#  variables from an input MC file. The input 
#  should be from the output of OscNtupleBuilder 
#  and contain generalized type, momentum, and dir
#  information
#


def SkipType (  mode , type ) :

    skip = 0 

    # this is done to maintain the 
    # fortran indexing withing OscTypes
    cmp =  (itype[0]-1) % g.nTypes
    cmp += 1 

    # skip upmu thru shower and upthro non-shower
    if ( \
         cmp == g.UpThruNonShower_mu  or \
         cmp == g.UpThruShower_mu     ) :  skip = 1 

    # skip neutral current events 
    if( abs(mode) > 30 ): skip = 1 

    return skip


def BuildNeighbors( FullNeighbors, nNeighbors = 20  ):
   #print 'Build Neighborts'
  
   dict = {}

   length = len( FullNeighbors )
   
   j = 0 
   for x in FullNeighbors :
   
     #default averaging
     nAve = nNeighbors 

     #number +/- in array to draw from
     # beginning
     nPM = ( j )*2  
     if( nAve > nPM  ) : nAve = nPM         
     if( nAve <  0   ) : nAve = 0          
     
     #number +/- in array to draw from
     # end 
     nPM = ( length - j - 1 ) * 2
     if( nAve > nPM  ) : nAve = nPM         
 
     neighbors = []
     for i in range ( -nAve/2 , nAve/2 + 1 , 1 ): 
        if( i == 0 ) : continue 
        neighbors += [ FullNeighbors[j+i].pnu ]

#    print j , neighbors 
     dict [ x.event ]  =  neighbors 
     #print 'event=',x.event
     #print 'neighbors=', neighbors

     j += 1 
 

   return dict 




# entry point for command line execution
if __name__ == "__main__":
   import commands
   import sys

   if ( len( sys.argv ) != 4 ): 
      print "Usage: " , sys.argv[0] , "[infile.root] [tree] [outfile.root]"
      sys.exit()

   infile     = str (sys.argv[1])
   treename   = str (sys.argv[2])
   outfile    = str (sys.argv[3])

   inf   = TFile( infile )
   ltree = inf.Get( treename )
   nEntries = ltree.GetEntries()
   print nEntries

   itype = array.array('i' ,    [0] )
   mode  = array.array('i' ,    [0] )
   ipnu  = array.array('i' ,    [0] )
   amom  = array.array('f' ,    [0] )
   dir   = array.array('f' ,  3*[0] )
   pnu   = array.array('f' ,    [0] )


   ltree.SetBranchAddress( 'itype' , itype ) 
   ltree.SetBranchAddress( 'ipnu'  , ipnu  ) 
   ltree.SetBranchAddress( 'pnu'   , pnu   ) 
   ltree.SetBranchAddress( 'dir'   , dir   ) 
   ltree.SetBranchAddress( 'amom'  , amom  ) 
   ltree.SetBranchAddress( 'mode'  , mode  ) 

   g = OscTypes
# rvw
   #nEntries = 50000 

   objects = [] 
   for j in range( 0, nEntries ): 
       ltree.GetEntry(j)

       if( j % 100000 == 0 ):  print "Entry Number: ", j

       if( SkipType( mode[0] , itype[0] ) == 1 ) : continue 

       zbin = g.GetZBin( -dir[2] , itype[0] )

       ## ( self , type, cz, amom, pnu, ipnu , ievent ):
       objects += [  EAveContainer(                \
                                     itype[0]    , \
                                     zbin        , \
                                     amom [0]    , \
                                     pnu  [0]    , \
                                     ipnu [0]    , \
                                     j
                                  )] 
   # --------- loop on tree entries has ended

   print "now sorting"
   result = sorted( objects, key=attrgetter('type' , 'zbin' , 'ipnu' , 'amom'  ) ) 

   min = int( ltree.GetMinimum('itype') )
   max = int( ltree.GetMaximum('itype') )

   print 'TYPES: ' , min, max

   RefList = [] 
   for type in range( min, max + 1 ):
     for z in range(0,10):
      for ip in ( -16, -14, -12 , 12 , 14 , 16 ) :
       subset = [ x for x in result if x.type == type and x.zbin == z and x.ipnu == ip ]
 
#      print type, z , ip , subset 
       
       # incase there is nothing 
       if( len( subset ) <= 0 ) : 
       #  print "   no subset info for  zbin: ", z , " type: ", type , " ip: ", ip 
          continue 

       print " zbin: ", z , " type: ", type , " ip: ", ip, subset[0] 

       dict = BuildNeighbors( subset ) 

       RefList += [ RefContainer( type, z , ip, dict ) ]
       

   # ----- loop on type, z, ip has ended 

   # Prepare output
   ouf   = TFile ( outfile      , 'recreate' )
   otree = TTree ( 'eaveFriend' , 'friend tree for eneryg averaging' ) 

   nlist  = array.array('i' ,      [0] )
   elist  = array.array('f' ,   20*[0] )
   dlist  = array.array('f' ,   20*[0] )
   mean   = array.array('f' ,      [0] )
   rms    = array.array('f' ,      [0] )
   parent = array.array('f' ,      [0] )

   # for direct comparison with osc3d's nonsense
   nlisthax = array.array('i' ,      [0] )
   elisthax = array.array('f' ,   20*[0] )
   rmshax   = array.array('f' ,      [0] )

   otree.Branch( 'nEAve'     ,  nlist  ,  'nEAve/I' )
   otree.Branch( 'NeighborE' ,  elist  ,  'NeighborE[nEAve]/F' )
   otree.Branch( 'DiffList'  ,  dlist  ,  'DiffList[nEAve]/F' )
   otree.Branch( 'parent'    ,  parent ,  'parent/F' )
   otree.Branch( 'Emean'     ,  mean   ,  'Emean/F' )
   otree.Branch( 'Erms'      ,  rms    ,  'Erms/F' )

   otree.Branch( 'nEAveHax'     ,  nlisthax  ,  'nEAveHax/I' )
   otree.Branch( 'NeighborEHax' ,  elisthax  ,  'NeighborEHax[nEAveHax]/F' )
   otree.Branch( 'ErmsHax'      ,  rmshax    ,  'ErmsHax/F' )

   print "Now filling output tree ... "
   for j in range( 0, nEntries ): 
       ltree.GetEntry(j)

       parent  [0] = pnu[0] 
       mean    [0] = 0.0
       rms     [0] = 0.0
       rmshax  [0] = 0.0
       nlist   [0] = 0
       nlisthax[0] = 0
       for i in range( 0, 20 ):
          elist   [i] = 0. 
          dlist   [i] = 0. 
          elisthax[i] = 0. 

       # skip neutral current event and upmu through events 
       if( SkipType( mode[0] , itype[0] ) == 1 ) :
          otree.Fill()
          continue


       zbin = g.GetZBin( -dir[2] , itype[0] )
       Ref = None 
       for x in RefList :
          if( x.type == itype[0] and x.zbin == zbin and x.ipnu == ipnu[0] ):
             Ref = x 
             break 
      
       if Ref is None : 
         print ' ----- Ref is None ----- '
         # When there are no neighbors to average 
         # foce nEAve  
         # from barfing 
         nlist[0]    = 1 
         nlisthax[0] = 1 
         otree.Fill()
         continue 
           
       nlist[0] = len( Ref.dict[j] )
       if( nlist[0] == 0 ):
          print "     CC event no neghbors: ", j , itype[0], zbin , ipnu[0], \
                     Ref.zbin, Ref.type , Ref.dict[j] , len( Ref.dict[j] ) 
 
       if( nlist[0] <= 0 ):
         # fill with ones to prevent ROOT 
         # from barfing 
         nlist[0]    = 1 
         nlisthax[0] = 1 
         otree.Fill()
         continue

       for i in range( 0, nlist[0] ) :
          elist[i] = Ref.dict[j][i]
          dlist[i] = pnu[0] - Ref.dict[j][i]
          mean[0] += elist[i] / float( nlist[0] ) 

       for i in range( 0, nlist[0] ) :
          rms[0] +=  ( elist[i] - pnu[0] )*( elist[i] - pnu[0] )
       rms[0] = sqrt( rms[0] / float(nlist[0]) ) 

 
       # compute haxxed variables to compare with 
       # osc3d scheme (ie only accept energies that with 0.5 or 2.0* pnu)
       q = 0 
       for i in range( 0, nlist[0] ) :
          if( Ref.dict[j][i] < 2.0*pnu[0] and Ref.dict[j][i] > 0.5*pnu[0] ) :
             elisthax[q] = Ref.dict[j][i]
             q += 1 
             nlisthax[0] += 1 

       rmshax[0] = 0.0
       for i in range( 0, nlisthax[0] ) :
          rmshax[0] +=  ( elisthax[i] - pnu[0] )*( elisthax[i] - pnu[0] )

       if ( nlisthax[0] > 0 ):
          rmshax[0] = sqrt( rmshax[0] / float(nlisthax[0]) ) 

       print 'rmshax=',rmshax[0]

       otree.Fill()
   # -------- looop on tree entries has ended 

    
   # finish up
   otree.Write("", TObject.kOverwrite)
   ouf.Close()

   inf.Close()


#
#  End of routine
#
############

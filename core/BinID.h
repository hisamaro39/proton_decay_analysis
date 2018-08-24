#ifndef _BinID_
#define _BinID_

#include <iostream>
#include <string>
// Simple base class for labelling bins in the analysis. 
//   A Tag represents the number of the bin as it appears
//   in the list of bins maintained by the profiler. 
//   Type is specified by the user, also an integer. A child
//   of this class should be made in the event.
//   there are external parameters controlling the bin labelling,
//   such as Energy, Momentum, Zenith Angle, etc.
//

class BinID
{
 public:
  BinID(){ 
    tag  = -1; 
    type = -1;
    _pEdges  = 0;
    _pBinary = 0;    
  }

  BinID( int d, int a , int b , const char * llabel = "none" ) 
  { 
       NEdges   = d; tag = a; type = b ; 
       _pEdges  = new double[NEdges];
       _pBinary = new int[NEdges];     
       Label    = llabel; 

  }

  ~BinID()
  { 
    if( _pEdges  != 0 ) delete [] _pEdges; 
    if( _pBinary != 0 ) delete [] _pBinary;
  }
  


  void SetEdge( int d, double x ) { _pEdges[d]  = x; }
  void SetDimN( int d, int N   ) { _pBinary[d] = N; }
  int  GetnEdges() { return NEdges ; }

  int    Tag()     { return tag; };     
  // tag specifies which bin number of the analysis the object represents //
  int    Type()    { return type; };

  std::string  GetLabel() { return Label; }

 
  double Edge( int d ) { return  _pEdges[d]; };
  int DimBinNum( int d ) { return _pBinary[d];};
 
  void Print(){ };
  
 protected:

  
  double *_pEdges;
  int  *_pBinary;     // Keeps track of which number of which dimension the bin is...

  int NEdges;
  int tag;      // Specify which bin this id belongs to //
  int type;     // What type of Bin we are labeling //
  
  std::string Label;  // What type of bin is this, for humans
  
};
#endif



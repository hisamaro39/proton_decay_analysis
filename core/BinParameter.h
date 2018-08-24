#ifndef _BinParameter_
#define _BinParameter_
#include <string>

// This class essentially specifies the spacing of Bins, and is mostly used by the profiler.

class BinParameter {
 public:

  BinParameter( int lType, int *BinsPerDim , int d , const char * llabel) { 
             Type = lType ; 
             nDimensions = d; 
             Label = llabel;
             __pEdges = new double*[nDimensions]; 
             _BinsPerDim = new int[nDimensions];
             for ( int i = 0 ; i < nDimensions ; i++ )
                 _BinsPerDim[i] = BinsPerDim[i];
     
             nBins = 1;
             for ( int i = 0 ; i < nDimensions ; i++ )
                 nBins *= _BinsPerDim[i];

             EdgeLabels = new std::string[nDimensions];
          };

  
  void SetOffset( int c ) { X = c;};


  void SetEdges(int d, double *p ) { 
             if (d < nDimensions )
             {
                int n = _BinsPerDim[d];
                __pEdges[d] = new double[ n ];
                for( int i = 0; i < _BinsPerDim[d] ; i++ )
                   __pEdges[d][i] = p[i]; 
             }
         };
  

  void SetEdgeLabel(int d, const char *llabel ) { 
             if (d < nDimensions )
             {
                EdgeLabels[d] = llabel ;
             }
         };


  double GetEdge( int d, int n ) { 
          if ( d < nDimensions ) 
             if ( n < _BinsPerDim[d] )
               return __pEdges[d][n];
               
          return -1.0;
  } //end getEdge



  int GetType()             { return Type; };
  int GetnBinsPerDim(int d) { if (d < nDimensions ) return _BinsPerDim[d]; return -1; };
  int GetnBins()            { return nBins; };
  int GetOffset()           { return X; };
  int GetnDims()            { return nDimensions; };
  
  std::string GetLabel()            { return Label; };
  std::string GetEdgeLabel(int d)   { if ( d < nDimensions ) return EdgeLabels[d] ; else return "none" ;}

 ~BinParameter() 
  {
     delete [] __pEdges;
     delete [] _BinsPerDim;
     delete [] EdgeLabels;
  }
     

 private:

  int nDimensions;      // Number of Bin Dimensions
  int nBins;
  std::string Label; 
  std::string * EdgeLabels; 

  int *_BinsPerDim; //collection of number of bins per dimension

  double **__pEdges; //edges per bin dimension 

  int Type;
  int X;      //This listing of Bin's Offset
     
};
#endif

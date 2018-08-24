#ifndef _profiler_
#define _profiler_

#include <vector>
#include <math.h>
#include <sstream>
#include <map>

#include "core/BinParameter.h"
#include "core/EventParser.h"
#include "core/BinID.h"
#include "core/ProfileSpace.h"
#include "core/GenericOscPoints.h"

#include "tools/CardReader.h"

#include "core/coreLibVer.h"

/*
  This class will serve as an Interface between "bin" objects 
  and the rest of the computation.
  so that they can be used for things 
  other than a 3-flavor oscillation calculation.
  
  profiler is desined to carry various seperate bits of the analysis 
  in one place. Noteably, the BinParameters, the BinID's.
  In this way the program (via SysMaestro) can interact with a commmon
  source to gain information about the Binning scheme. 
  Note that SysMaestro will typically contain
  the Errors, while the profiler contains the Bins 
  (or at least a list of them and their labels).
  
  Because the profiler has this information it is used also to "Bin Events", 
  that is taking an 
  EventParser argument, determine which bin of the analysis a particular event 
  will contribute.
	 20060315
*/

class BinID;

class profiler {

 public:
  profiler(int , bool k = false );
  profiler( CardReader * );
  virtual ~profiler();
  
  
  virtual BinID* GetBinId( int );		
  /* return bin_ID corresponding to master bin number mostly useful for printing statements
     and other forms of debugging.*/
  
  virtual int BinEvent( EventParser* E);
  /* Return the master bin number for a given (parsed) Event */
  
  int GetnBins() { return nBins; };
  /* Enable external classes access to the number of Bins in the analysis */
  
  BinID* associate_bin( int type ); /* assign a created id to a bin */
  BinID* associate_bin( ); /* assign a created id to a bin */
  void init_association(); /* start bin counters at zero */
  

  // Get bin corresponding to a list of paramenters 
  // and a corresponding event type
  int GetBin(std::vector<double>, int );

  // Get bin corresponding to a list of parameter offsets 
  // and a corresponding event type. Ie the 
  // bin at (type, n1, n2,...nX) 
  // where n_i is the n^th bin number of dimension i
  int GetBinLocale(std::vector<int>, int );

  std::vector<BinParameter*> GetBinParameters() { return _vpBP; };
  /* return the vector of pointers to BinParameters 
     -- useful mainly for printing in situation where there is more than one type of event	
     BE CAREFUL this makes _vpBP basically public! This is bad form and will be corrected in
     future versions of the code									*/
 
  int GetnPoints() { return nPoints; } ; 
  int GetBinType(int);
  std::string GetBinLabel(int);

  void SetStartPoint(int a ){  StartPoint = a; }
  void SetStopPoint(int a  ){  StopPoint = a;  }

  int GetStartPoint( ){  return StartPoint; }
  int GetStopPoint( ){   return StopPoint;  }

  virtual int GetBinLocale( int a, int b, int c ){ return -1;}
  void TryToCreatePoints(){  if( !kTenCreated ) CreatePoints(); }

  void GenerateGenericBins( int a );

  // moved here to allow for simple fitting 
  void AddPointType( double , double , int , int , const char * name = NULL, int wrappable = 0, double * list = NULL);

  ProfileSpace& GetProfilePoint(int i) { return *_vpPS[i]; }

  std::vector< GenericOscPoints * > GetOscPointTypes() { return  OscPointTypes; }

  unsigned GetnSensPoints() { return _vSensitivityPoints.size() ; }  
  int      GetSensPoint( int i ) { return _vSensitivityPoints[i] ; }

  int GetStepsBetweenPoints(int,int);
  std::vector<int> GetProfilePointNeighbors(int i);
  CardReader * GetCardReader() { return cr ; }

  void LoadBinsFromCard( CardReader * lReader );

  bool GetAllowNegatives()  { return kAllowNegativeBins; }

 protected:
  CardReader   * cr ;

  /* the paramaters defining the binning scheme */
  ProfileSpace * SpaceMother;

  std::map< int, ProfileSpace* >   _vpPS;
  std::vector<BinID*>   _vpID;		/* vector for bin id informaition */
  std::vector<BinParameter*> _vpBP;		
 
  int nPoints;		/* Number of Different Oscillation Points */
 
  BinID	*_dummyID;	/* ID to return in case of error */
  int nBins;		/* total number of bins */
  int nAssociated;	/* Count of ID's which have been assigned(associated) to a bin */

  int NBinTypes;

  //This is the BinLayout Portion of the Profiler 
  void CreateBins();		
  void AddBinType( int TypeID, int NDimensions , int *BinsPerDimension , TString Label = "none" );
  void AddEdges( int TypeID, int Dim, double * , TString Label = "none" );

  void InitPoints(){ OscPointTypes.clear(); }
  void CreatePoints();

  std::vector< int > _vSensitivityPoints;

  std::map<int, BinParameter*>	_mpBP;
  std::map<int, int>		_mnBins;
  std::map<int, int >           _nA;		/* array to count which bin types have already been created */
  
  std::vector< GenericOscPoints* > OscPointTypes; 
 
  std::map < std::string , double * > BinEdges ;
  int StartPoint;
  int StopPoint;
  bool kTenCreated;      // used to prevent points from being made more than once
  bool kSensPointsExist; // Should other BinContainers be loaded to handle
                         // sensitivity points of Feldman Cousins points that 
                         // exist outide of the range of points covered by 
                         // this instance a program 
  bool kVerbose;
  bool kAllowNegativeBins; 

  coreLibVer * clvp;
};



#endif


#include "core/profiler.h"
#include <iostream>
#include <algorithm>
using namespace std;
profiler::profiler(int a , bool k)
{
  kVerbose = false;
  kAllowNegativeBins = k ;

  // Initialize the profiling vectors //
  NBinTypes = 0;
  _dummyID = new BinID() ;// defaults to dummy values type = -1!
  nBins = 0;
  nPoints = 1;

  // These are used for parallel computing, and MUST be set before
  // a call to CreatePoints() is made
  StartPoint = 0;
  StopPoint = nPoints;

  kTenCreated = false;

  coreLibVer * clvp = coreLibVer::GetLibrary();
  LibraryVerMaster::GetMaster()->AddLibrary( clvp ); 
}


profiler::profiler( CardReader * lReader )
{
  cr = lReader;
  nBins = 0;
  NBinTypes = 0;
  _dummyID = new BinID() ;// defaults to dummy values type = -1!

  // These are used for parallel computing, and MUST be set before
  // a call to CreatePoints() is made
  nPoints = 1;
  StartPoint = 0;
  StopPoint = nPoints;

  kVerbose = lReader->FindKey("profiler_verbose" );
  cout << "kVerbose = " << kVerbose << endl;

  kAllowNegativeBins = lReader->FindKey("profiler_negative_bins" );

  double Start;
  double Stop;
  int    kIsLog;
  int    nSteps;   
  int    nEntries;
  int    kWrappable;
  std::map<std::string, std::string > parms;
  std::map<std::string, std::string >::iterator _i;
  lReader->BuildListOfStrings( "parm_", parms );

  double *val;
  std::string::size_type loc;
  std::string::size_type loc2;
  TString  name;
  std::string Value;

  for( _i = parms.begin() ; _i != parms.end() ; _i++ )
  {
    
     loc  = _i->first.find("parm_", 0 , 5);
     if( loc == std::string::npos ) continue; 
     
     //extract everything after the suffix
     name  = _i->first;
     name.ReplaceAll("parm_", 5,"",0 );
     Value = _i->second;
     nEntries = lReader->ParseArray( Value , "," , &val );

     Start  = val[1];
     Stop   = val[2];
     kIsLog = (int) val[3]; 
     nSteps = (int) val[4]; 

     if( nEntries > 4 ) kWrappable = (int) val[5];
     else               kWrappable = 0;

     // IsLog = 2 is for hand specified parameter points
     // listed after the wrappable parameter
     if ( kIsLog != 2 ) AddPointType( Start, Stop, kIsLog, nSteps, name.Data(), kWrappable ) ;
     else               AddPointType( Start, Stop, kIsLog, nSteps, name.Data(), kWrappable, &val[6] ) ;
     
     delete [] val;

  }

 
  ///////////////////////
  //
  //  Store information for sensitivity points
  //
  ///////
  kSensPointsExist = false;
  std::string plane;
  bool kSens = lReader->GetKey("SensitivityArray" , plane     );
  
  if( kSens )
  {
     double * lArray;
     lReader->ParseArray( plane  ,  "," ,  &lArray  );
     int nSensitivity = (int) (lArray[0]-1);

     for( int i = 1 ; i < nSensitivity ; i++ )
     {
        _vSensitivityPoints.push_back( (int) lArray[i] );
        std::cout << "profiler::profile will store MC additionally at point: " 
                  <<  (int) lArray[i] << std::endl;
        if( i == 1 )
          kSensPointsExist = true;
     }
     delete [] lArray;
  }


// The Points will be created later
  kTenCreated = false;

  coreLibVer * clvp = coreLibVer::GetLibrary();
  LibraryVerMaster::GetMaster()->AddLibrary( clvp ); 

  return;
}

void profiler::LoadBinsFromCard( CardReader * lReader )
{

  std::map<std::string, std::string > ListOfTypes;
  std::map<std::string, std::string > ListOfAxisLabels;
  std::map< int       , std::string > OrderedListOfTypes;
  std::map<std::string, std::string > ListOfAxes;
  std::map<std::string, std::string >::iterator _i;
  std::map< int       , std::string >::iterator _oi;
  std::map<std::string, std::string >::iterator _j;
  TString sType;
  TString sAxisNumber;
  TString sAxisLabel;

  unsigned TypeNumber = 0;
  unsigned AxisNumber = 0;
  unsigned nAxes      = 0;
  unsigned size       = 0;
  int      BinsPerDim[10];
  double  * val;
  std::string Value;
  std::string AxisTag[ 10 ];
  std::string AxisLabel[ 10 ];

  lReader->BuildListOfStrings( "BinningType_", ListOfTypes       );

  

  std::string::size_type loc;
  std::stringstream ss;
  for( _i = ListOfTypes.begin() ;  _i != ListOfTypes.end() ; _i++ )
    OrderedListOfTypes[ (int) atoi( _i->second.c_str() ) ] = _i->first; 

  for( _oi = OrderedListOfTypes.begin() ;  _oi != OrderedListOfTypes.end() ; _oi++ )
  {
     
     loc  = _oi->second.find("BinningType_", 0 , 12);
     if( loc == std::string::npos ) continue; 
     
     //extract everything after the suffix
     sType = _oi->second;
     sType.ReplaceAll("BinningType_", 12,"",0 );

     //Store the Number we will globally assign to this type
     TypeNumber = (unsigned) ( _oi->first ); 

     // Now we will build a list of strings containing all of the 
     // axes for our new type. Then we can register the bins 
     ss.str(""); ss << "BinAxis_" << sType.Data() ;
     if (kVerbose) 
         std::cout << "Building Bins for : " << sType.Data() << " , " << ss.str() << " " << TypeNumber << std::endl; 
  
     nAxes = 0;
     lReader->BuildListOfStrings( ss.str().c_str() , ListOfAxes );

     // Grab the Axislabels (if they exist )
     ss.str(""); ss << "AxisLabel_" << sType.Data() ;
     lReader->BuildListOfStrings( ss.str().c_str() , ListOfAxisLabels  );

     for( _j = ListOfAxes.begin() ;  _j != ListOfAxes.end() ; _j++ )
     {
   
        sAxisNumber = _j->first;
        // expect the last "_" to preceed an integer specifying the axis number
        loc = sAxisNumber.Last('_');  
        sAxisNumber.Remove( 0, loc+1 );
        AxisNumber = atoi( sAxisNumber.Data() );

        if ( AxisNumber < 0 || AxisNumber > 10 ) 
        {
           std::cout << " profiler::LoadBinsFromCard Error - bad axis number : " << AxisNumber << std::endl;
           std::cout << "   BinAxis : " << _j->first << " will be skipped " << std::endl;
           std::cout << std::endl;
           continue ;
        }
        nAxes ++ ;

        Value = _j->second;
        int nEntries = lReader->ParseArray( Value , "," , &val );
        if (kVerbose)
            std::cout << " nEntries: " << nEntries << " " << _j->first ; 

        // one less bin than entries due to final edge
        BinsPerDim[ AxisNumber ] = nEntries - 1   ;
        AxisTag   [ AxisNumber ] = _j->first      ;

        ss.str(""); ss << "_" << AxisNumber;
        AxisLabel[ AxisNumber ] = ss.str().c_str() ;
        // Find the Correct Label for the axis, if it exists
        if( ! ListOfAxisLabels.empty() )
        {
            ss.str(""); ss << "AxisLabel_" << sType.Data() << "_" << AxisNumber;
            if( ListOfAxisLabels.count( ss.str().c_str() ) )
                AxisLabel[ AxisNumber ] = ListOfAxisLabels[ ss.str().c_str() ];
        } 

        // remove the antiquted padding at the end of the array 
        size  = (unsigned) val[0] - 1;

        // adjust the array and save in memory
        BinEdges[ _j->first ] = new double [ nEntries ] ;
        if (kVerbose) std::cout <<  "    " ;
        for( unsigned k = 1 ; k < size ; k++ )
        {
           BinEdges[ _j->first ][ k-1 ] = (double) val[k] ; 
           if (kVerbose) std::cout <<  " " << val[k] ;
        }
        if (kVerbose) std::cout << std::endl;

     } // end of loop on BinAxes


     std::string Label;
     if( kVerbose )
     {
        std::cout << " Registering Bins for " << _oi->second << " Type: " << TypeNumber << std::endl;
        std::cout << "   using  : "  << nAxes << " axes. " << std::endl;
     }

     // Now register the bins
     AddBinType( TypeNumber  , nAxes , BinsPerDim, sType );
     int nBinsAdded = 1;
     for (unsigned q = 0 ; q < nAxes ; q++ )
     {
        Label = AxisTag[ q ];
        nBinsAdded *= BinsPerDim[q] ;
        AddEdges( TypeNumber , q , BinEdges[Label] , AxisLabel[q] ); 
        if( kVerbose )
          std::cout << "  Adding axis " << q << " with " << BinsPerDim[q] <<  " Label: " << AxisLabel[q] <<std::endl; 
     }
     if( kVerbose )
        std::cout << " Added " << nBinsAdded << " for Type " << sType.Data() << std::endl;  

     // separate storage provided by "BinEdges" delete the 
     // temporary storage
     delete [] val;

   }// end of list of bin types
     


   CreateBins();
   std::cout <<"profiler::LoadBinsFromCard Created " << nBins << " bins from the input card" << std::endl; 
   
}






BinID* profiler::GetBinId( int i )
{
  return _vpID[i];
}

void profiler::init_association()
{
  nAssociated = 0;
  std::map<int,int>::iterator i;
  if( kVerbose )
     std::cout << "profiler::profiler NBinTypes "  << NBinTypes << std::endl;

  for ( i= _nA.begin(); i != _nA.end() ; i++ )
    i->second = 0;
}

BinID* profiler::associate_bin( int type )
{
  int i=0;

  // this is an error! we have created too many total bins! //
  if ( nAssociated > nBins ) return _dummyID;       

  // this is an error! we have created more bins of a given type than we planned for! //
  if ( _nA[type] > _mnBins[type] ) return _dummyID; 
  
  int vpIDsize = _vpID.size();
  while( _vpID[i]->Type() != type &&  i < vpIDsize ) 
    i++;
  
  nAssociated++;

  // _nA[type] records the current number of bin ID's created for a given type //
  _nA[type]++;
  
  return _vpID[i+_nA[type]-1];
}

BinID* profiler::associate_bin()
{
   if ( nAssociated >= nBins ) return _dummyID;  

   // this is an error! we have created more bins of a given type than we planned for! //
   if ( nAssociated > 0 )
     if ( _nA[ _vpID[nAssociated-1]->Type()] > _mnBins[ _vpID[nAssociated-1]->Type() ] ) 
        return _dummyID; 
        
   if ( nAssociated > 0 )
      _nA[ _vpID[nAssociated-1]->Type() ]++;

   nAssociated++;

   return _vpID[ nAssociated-1];
}


int profiler::GetBinType(int n )
{
   return _vpID[n]->Type();
}




void profiler::AddBinType( int TypeID, int NDimensions, int *BinsPerDimension , TString Label  )
{

   int i;
   int Bins = 1;

   TString PassLabel;
   if ( Label != "none" ) PassLabel = Label;
   else 
   {
     std::stringstream ss;
     ss.str(""); ss << "BinType" << TypeID ;
     PassLabel = ss.str().c_str();
   }
   
   _mpBP[ TypeID ] = new BinParameter( TypeID, BinsPerDimension, NDimensions, PassLabel.Data() );
   _vpBP.push_back( _mpBP[ TypeID ] );

   for( i = 0 ; i < NDimensions; i++ )
      Bins *= BinsPerDimension[i];

   if ( _mnBins.count( TypeID ) == 0 )
   {   
      NBinTypes++;
      _nA[ TypeID ] = 0;
   }


   _mnBins[ TypeID ] = Bins;


  if( kVerbose )
    std::cout << "profiler::AddBinType " << Bins << " " << TypeID << std::endl;   
}



void profiler::AddEdges( int TypeID, int Dim, double *Edges, TString Label  ){

   _mpBP[TypeID]->SetEdges    ( Dim, Edges );
   _mpBP[TypeID]->SetEdgeLabel( Dim, Label.Data() );
}


// 
//    CreateBins():
//    This may seem like a very odd fucntion to one who hasn't seen it before. The main
// idea is to createa sequential listing of binID's based upon the users input to the derived
// contsructor of _Name_profiler.cc. The user may add bins in any order, so long as they have a unique
// type, specified in global.h. This function creates the order of the bin, based on the order specified
// in global.h.  
//
//    A call to AddBinType enables one to add a bin which is indexed by an arbitrary number of dimensions.
//  
//   Examples: 
//      - the SK ATMPD analysis bins events in  ( Lepton Momentum , Zenith Angle ) , ie 2 dimensions
//      - the T2K analysis bins events in  ( Energy ), ie just 1 dimension;
//
//   Because the number of dimensions is decided by the user, this routine has to decide automatically how to create
//the binID's (which in turn create the bins themselves ). The Convention is the last dimension changes the fastest, the
//first dimension changes the slowest.
//
//
//   How it works.
//
//   Suppose you have a BinType with 4 dimensions, with ( 4, 2, 7, 6 ) bins in each dimension, for a total of 336 Bins
//    for this type. To create the 200th bin:
//
//       200 = (2*7*6)w + (7*6)x + 6y + z + 1    , here (w,x,y,z) represent which number of which dimension this
//                        bin represents.
//
//      Please note the factor multiplying each of w,x,y,z!
//
//      Starting with the slowest index, w =  [ (200-1)/(2*7*6) ] , where [A] is the greatest integer <= A.
//          
//            then x =  [  ((200-1) % (2*7*6)) / (7*6) ], and so on until all (w,x,y,z) 
//                               have been specified.
//
//   All the ID's are created by looping over the total number of bins for a given type and following this algorithm.
//


void profiler::CreateBins()
{

   int BinCount = 0;
   int Type;

   std::map<int, BinParameter*>::iterator __i;

   BinID *_tmp;

   int *DimFactors;
   int *Dims;
   int *BinBinary;
   double Edge;
   int  mynBins;
   std::string Label;


   int Quotient, Remainder, Hack;
   int i, j,  NDims, BinNumber;

   // Loop over BinParameters
   // The map automatically sorts the bins by key, ie, the order in global.h!
   for ( __i = _mpBP.begin() ; __i != _mpBP.end() ; __i++ )
   {
      //set How Far along in the list of bins is this type
      __i->second->SetOffset(BinCount);

      //get number of Dimensions from the BinParameter
      Type  = __i->second->GetType();
      fprintf(stderr,"%i /n",Type);
      NDims = __i->second->GetnDims(); 
      Label = __i->second->GetLabel();  

      DimFactors  = new int[NDims];
      Dims        = new int[NDims];
      BinBinary   = new int[NDims];

      for( i = 0 ; i < NDims ; i ++ )
         Dims[i] = __i->second->GetnBinsPerDim(i);   // Get Number of Bins per Dimension;

   
      for( i = 0 ; i < NDims ; i ++ )
      {
         DimFactors[i] = 1;
         if ( i < NDims -1 )
            for( j = i + 1; j < NDims ; j++ )
               DimFactors[i] *= Dims[j];
      }         

      mynBins = _mnBins[Type];

      // BinNumber here refers to the bin Number *within* the type of bin
      // eg a BinType with 2 dimensions, with (6,8) divisions respectively
      // will have BinNumbers = (1,48) to choose from.
      // this is a small break from normal indexing
      for ( BinNumber = 1 ; BinNumber <= mynBins ; BinNumber++ )
      {
         for( i = 0 ; i < NDims  ; i++ )
         {

            if ( i == 0 )
              Hack = BinNumber - 1;              // the indexing adjustment is made here

            Quotient = int( floor( double( Hack) / double(DimFactors[i]) ) );
            Remainder = Hack % DimFactors[i];   
      
            BinBinary[i] = Quotient;

            Hack = Remainder;
         }

         _tmp = new BinID( NDims, BinCount , Type , Label.c_str() );
         for( i = 0 ; i < NDims ; i++ )
         {
            Edge = __i->second->GetEdge(i, BinBinary[i] );   
            _tmp->SetEdge( i, Edge );
            _tmp->SetDimN( i, BinBinary[i]);
         }

         _vpID.push_back( _tmp );

         BinCount++;

      } // End Loop over BinNumbers

      delete[] DimFactors;
      delete[] Dims;
      delete[] BinBinary;

   } // End of Loop over BinTypes

   nBins = BinCount;
}

void profiler::AddPointType(double Start, double Stop, int isLog, int nSamp, const char * name , int kWrappable, double * list )
{

    nPoints *= nSamp;
    std::stringstream ss;
    
    if ( name != NULL ) 
      ss << name ;
    else
       // if no name is specified use the parameter number
      ss << OscPointTypes.size() ;
       
    if ( list != NULL )
    {
       GenericOscPoints * ptr = new GenericOscPoints( nSamp, ss.str().c_str() , kWrappable );
       for( unsigned i = 0 ; i < nSamp ; i++ )
         ptr->SetPoint( i , list[i] );

       OscPointTypes.push_back( ptr );
    }
    else 
       OscPointTypes.push_back( new GenericOscPoints( Start, Stop, isLog, nSamp, ss.str().c_str(), kWrappable ) );
    
}





// This function creates a list of points to test in the model space of arbitrary dimension
// the scheme is identical to CreateBins above
// 
//         -rvw

void profiler::CreatePoints()
{

  int OscSize = (int)OscPointTypes.size();
  int * Coefficients = new int[ OscSize ];
  int * Index = new int[ OscSize ];
  int Remainder;
  double * PointValues = new double [OscSize]; 
  double localValue;
  double lStart, lStop, lnPoints; //l is for local
  int i,j;
  int pointcounter = 0;

  for ( j = 0 ; j < OscSize; j++ ){
    Coefficients[j] = 1;
    for ( i = j+1 ; i < OscSize; i++ )
      Coefficients[j] *= OscPointTypes[i]->nSamples;
  }    

   // register the variable names to enable calls to ProfileSpace::Get("blah");
   SpaceMother = new ProfileSpace( -100 , OscSize , new double[OscSize] );
   for( j = 0 ; j < OscSize ; j++ )
     SpaceMother->Register( OscPointTypes[j]->Name.c_str() , j , kVerbose ); 
    

// Create the list of points in profiling space: //
  for ( i = StartPoint ; i < StopPoint ; i++ )
  {

     //determine the idices of this point...
     Remainder = i;
     for( j = 0 ; j < OscSize ; j++ )
     { 
       Index[j] = (int)  floor( (double) Remainder / (double) Coefficients[j] )  ;
       Remainder = i % Coefficients[j]; 
     }

     // no these loops do not need to be separate, but are done so for clarity only
     for( j = 0 ; j < OscSize ; j++ )
     {
        lStart   = OscPointTypes[j]->Start;
        lStop    = OscPointTypes[j]->Stop;
        lnPoints = (double) OscPointTypes[j]->nSamples;
 
        localValue =  OscPointTypes[j]->GetPoint( Index[j] ); 
        PointValues[j] = localValue; 
     }
     _vpPS[i] =  new ProfileSpace( i, OscSize, PointValues ) ;

     pointcounter++;
  }// end of loop on points...


  if( kVerbose )
     std::cout << "profiler::CreatePoints " 
               << pointcounter << " Analysis Profile points have been Created " << std::endl;

  kTenCreated = true;

  delete [] Coefficients; 
  delete [] Index; 
  delete [] PointValues;

}

int profiler::GetStepsBetweenPoints(int a,int b)
{
  if(a>StopPoint or b>StopPoint or a<StartPoint or b<StartPoint)
    return -1;

  int diff=abs(a-b);

  std::vector<std::pair<int,int> > stepsVec;
  if(OscPointTypes.back()->kWrappable>0)
    {
      stepsVec.push_back(std::make_pair(OscPointTypes.back()->nSamples,OscPointTypes.back()->nSamples));
      }
  else
    {
      stepsVec.push_back(std::make_pair(OscPointTypes.back()->nSamples,-1));
    }
  for(int iOscParm=(int)OscPointTypes.size()-2;iOscParm>=0;iOscParm--)
    {
      if(OscPointTypes[iOscParm]->kWrappable>0)
	{
	  stepsVec.push_back(std::make_pair(stepsVec.back().first*OscPointTypes[iOscParm]->nSamples,OscPointTypes[iOscParm]->nSamples));
	}
      else
	{
	  stepsVec.push_back(std::make_pair(stepsVec.back().first*OscPointTypes[iOscParm]->nSamples,-1));
	}
    }
  
  int steps=0;
  std::vector<std::pair<int,int> >::iterator it;
  for(it=stepsVec.end()-1;it!=stepsVec.begin()-1;it--)
    {
      int theseSteps=diff/(it->first);
      if(it->second>0)
	{
	  theseSteps=std::min(theseSteps,it->second-theseSteps);
	}
      steps+=theseSteps;

      diff=diff%(it->first);

    }
  steps+=diff;


  return steps;

}

std::vector<int> profiler::GetProfilePointNeighbors(int i)
{

  std::vector<int> ret;
  int step=1;
  for(int iOscParm=(int)OscPointTypes.size()-1;iOscParm>=0;iOscParm--)
    {
      int ptu=i+step;
      int ptd=i-step;
      step*=OscPointTypes[iOscParm]->nSamples;
      //      std::cout<<" iOscParm: "<<iOscParm<<std::endl;
      //std::cout<<" i: "<<i<<" ptu: "<<ptu<<" ptd: "<<ptd<<" start: "<<StartPoint<<" stop: "<<StopPoint<<std::endl;
      if(ptu<StopPoint) 
      	{
	  //  std::cout<<"ptu val: "<<_vpPS[ptu]->Get(OscPointTypes[iOscParm]->Name.c_str())<<" i val: "<<_vpPS[i]->Get(OscPointTypes[iOscParm]->Name.c_str())<<std::endl;

	  //	  if(_vpPS[ptu]->Get(OscPointTypes[iOscParm]->Name.c_str())>_vpPS[i]->Get(OscPointTypes[iOscParm]->Name.c_str()))

	  if(_vpPS[ptu]->Get(OscPointTypes[iOscParm]->Name.c_str())>_vpPS[i]->Get(OscPointTypes[iOscParm]->Name.c_str()))
	  ret.push_back(ptu);
	  else if(OscPointTypes[iOscParm]->kWrappable>0)
	    {

	      ptu-=step;
	      ret.push_back(ptu);

	    }
	}
      else if(OscPointTypes[iOscParm]->kWrappable>0)
	{
	  ptu-=step;
	  ret.push_back(ptu);
	}
      if(ptd>=StartPoint)
	{
	  if(_vpPS[ptd]->Get(OscPointTypes[iOscParm]->Name.c_str())<_vpPS[i]->Get(OscPointTypes[iOscParm]->Name.c_str()))
	    ret.push_back(ptd);
	  else if(OscPointTypes[iOscParm]->kWrappable>0)
	    {

	      ptd+=step;
	      ret.push_back(ptu);

	    }


	}
      else if(OscPointTypes[iOscParm]->kWrappable>0)
	{
	  ptd+=step;
	  ret.push_back(ptd);

	}

      
      if(find(ret.begin(),ret.end()-1,ret.back())!=ret.end()-1) ret.pop_back();
    }

  return ret;

}

profiler::~profiler()
{
      
   delete _dummyID; 

   for ( int i = 0 ; i < _vpID.size() ; i++ )
      delete _vpID[i];
   
   std::map<int, ProfileSpace*>::iterator _i;
   for ( _i = _vpPS.begin() ; _i != _vpPS.end() ; _i++ )
      delete _i->second;


   for ( int i = 0 ; i < _vpBP.size() ; i++ )
      delete _vpBP[i];

// _vpBP and _mpBP hold the same objects
// map<int,BinParameter*>::iterator __i;
// for( __i=_mpBP.begin(); __i != _mpBP.end() ; __i++ )
//    delete    __i->second;   

   map<std::string, double *>::iterator ___i;
   for( ___i=BinEdges.begin(); ___i != BinEdges.end() ; ___i++ )
      delete [] ___i->second;   

   delete SpaceMother;
         
  for ( unsigned i = 0 ; i < OscPointTypes.size() ; i ++ )
    delete OscPointTypes[i];
}


void profiler::GenerateGenericBins( int nGen )
{
   int i,n;
   int nDims = 1;
   BinID * _tmp;
   int BinCount = 0;
   double BinEdges[ nGen ] ;

   _vpID.clear();
   for ( n = 0 ; n < nGen ; n++ )
   {
     _tmp = new BinID( nDims, n , n );
     BinEdges[n] = double(n); 

     for( i = 0 ; i < nDims ; i++ )   
     {
       _tmp->SetEdge( i , double(n) );
       _tmp->SetDimN( i , n     );
     }

     _vpID.push_back( _tmp );
   }
   
   int TypeID = 0 ;
   int BinsPerDimension [1] = { nGen };
   AddBinType( TypeID, nDims , BinsPerDimension , "GenericBins"  );
   AddEdges( TypeID , 0 , BinEdges ); 

   nBins = nGen;
}



int profiler::GetBin( std::vector<double> VarList, int type )
{
   if( ! _mpBP.count( type ) ) return -1;

   int nDimensions = _mpBP[type]->GetnDims();
   unsigned i,j,k;


   //for (unsigned q = 0 ; q < VarList.size() ; q++ )
   //  std::cout << " Var" << q << " " << VarList[q] << " " << type << std::endl;
   //std::cout << std::endl;   

   // We will allow the input list to be larger 
   // than the current bin type's allowed dimension
   if( VarList.size() < nDimensions )
   {
       std::cout << "Profiler::GetBin Error. Variable list is too small "  
                 << " size: " << VarList.size() << " expected " << nDimensions << endl 
                 << " for dimensions specified by bin type " << type << ". Exiting."
                 << std::endl ; 
       exit(-1);
   } 

   int DimFactors  [nDimensions];
   int LocalDimNum [nDimensions];
   for( i = 0 ; i < nDimensions-1; i++ )
   {
      DimFactors[i] = 1;
      for( j = i + 1; j < nDimensions ; j++ )
         DimFactors[i] *=  _mpBP[type]->GetnBinsPerDim( j );
   }
   DimFactors[nDimensions -1 ] = 1;


   int X  = _mpBP[type]->GetOffset();
   int SubBinOffSet;
   int BaseAddress;
   int ThisAddress;
   int NextAddress;

// Illustrative examples for bin addressing 
// in the linearalized schemne:
//
// If you have two dimensional bins
//  (Momentum p, and Zenith Angle Z )
//  with Np bins in p, and Nz bins in Z:
// 
//  The p bin edge ordering is:
//  p ->  X + p*Nz
//
//  The z bin edge ordering is:
//  z ->  X + p*Nz + z 

//
//  If three dimensional 
//  (Momentum p, Zenith Angle Z , Cats c )
//  p -> X + p*Nz*Nc
//  z -> X + p*Nz*Nc + z*Nc
//  c -> X + p*Nz*Nc + z*Nc + c
//
//  Below:
//  DimFactors for dimension1 (p) is Nz*Nc
//  DimFactors for dimension2 (z) is Nc
//  DimFactors for dimension3 (c) is 1
//

   int nBins;

   // loop over all of the bin dimensions
   for( i = 0 ; i < nDimensions ; i++ )
   {
      LocalDimNum[i] = -1;

      // compute the base offset for this bin types i_th dimension
      // Start at the global offset
      BaseAddress = X;
      if( i > 0 )
         for( k = 0 ; k < i-1 ; k ++ )
            BaseAddress += LocalDimNum[k]*DimFactors[k];        

      nBins = _mpBP[type]->GetnBinsPerDim( i );
      for( j = 0 ; j <  nBins ; j++ )
      { 

         // increment by one for the current bin edge
         ThisAddress = BaseAddress + j * DimFactors[i];        

         // increment by one more for the next bin edge
         NextAddress = ThisAddress + DimFactors[i];        

         // we are between two edges
         if( j != nBins -1 )
            if( VarList[i] >= _vpID[ ThisAddress ]->Edge(i) && VarList[i] < _vpID[ NextAddress ]->Edge(i) )  
            { 
               LocalDimNum[i] = j ; 
               break;
            }

         if( j == nBins -1 )
            // we are in the last bin 
            if( VarList[i] >= _vpID[ ThisAddress ]->Edge(i) )
            {
               LocalDimNum[i] = j ; 
               break;
            }

      }// end of loop on number of bins in this dimension


      // we reached the end of the bin edges 
      // and we didn't fall in any of the defined bins -->Error!
      if( LocalDimNum[i] < 0 ) 
         return -1;  
   

   } // end of loop on dimensions


   ThisAddress = X;
   for( k = 0 ; k < nDimensions ; k ++ )
      ThisAddress += LocalDimNum[k]*DimFactors[k];        

   return ThisAddress;

}

int profiler::GetBinLocale( std::vector<int> Indices , int type )
{

   int nDimensions = _mpBP[type]->GetnDims();
   unsigned i,j;

   // We will allow the input list to be larger 
   // than the current bin type's allowed dimension
   if( Indices.size() < nDimensions )
   {
       std::cout << "Profiler::GetBin Error. Indices list is too small " 
                 << " for dimensions specified by bin type " << type << ". Exiting."
                 << std::endl ; 
       exit(-1);
   } 

   int DimFactors  [nDimensions];
   int LocalDimNum [nDimensions];
   for( i = 0 ; i < nDimensions -1; i++ )
   {
      DimFactors[i] = 1;
      for( j = i + 1; j < nDimensions ; j++ )
         DimFactors[i] *=  _mpBP[type]->GetnBinsPerDim( j );
   }
   DimFactors[ nDimensions -1 ] = 1;

   int X  = _mpBP[type]->GetOffset();
   int ThisAddress = X;
   for( i = 0 ; i < nDimensions ; i++ )
   {
     // bounds checking on indices
     // fail if a requested bin does not exist
     if ( Indices[i] < _mpBP[type]->GetnBinsPerDim(i ) )
      ThisAddress += Indices[i] * DimFactors[i];
     else 
     {
      ThisAddress = -1;
      break;
     }
   }

   // check that we have not found a bin 
   // outside the range of this Type 
   if ( ThisAddress > 0 ) 
     if( type != GetBinId( ThisAddress )->Type() )
        return -1;

   return ThisAddress;
}

int profiler::BinEvent( EventParser* E )
{
   std::vector<double> vars;
   E->GetEventBinValues( vars );



   return GetBin( vars , E->GetType() );

}





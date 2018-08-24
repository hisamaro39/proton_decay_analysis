#ifndef _FileRecord_
#define _FileRecord_

#include <sstream>
#include "TString.h"
#include "TChain.h"
#include "tools/CardReader.h"

class FileRecord 
{

  public:

     // FilePath, Instance, nPoints, nFiles
     FileRecord( const char * , int , int , int ); 

     FileRecord( );

     TString BuildFileInstanceString (int p){return BuildFileInstanceString(FilePath,p);}
     TString BuildFileInstanceString   ( const char * , int  ); 
     TString BuildOutputInstanceString ( const char * s , int instance,  int nfiles  );

     int     GetInstanceFromPoint      ( int Point    , int npoints  = -1 , int nfiles  = -1  );
     int     GetPointOffsetFromInstance( int Point    , int instance = -1 , int npoints = -1  , int nfiles = -1 );
     void    GetSeedAndBlossom         ( int & Seed   , int & Blossom    , int instance = -1 , int npoints = -1 , int nfiles = -1 ) const;

     void    GetSeedAndBlossomOld         ( int & Seed   , int & Blossom    , int instance = -1 , int npoints = -1 , int nfiles = -1 ) const;
     void    GetSeedAndBlossomNew         ( int & Seed   , int & Blossom    , int instance = -1 , int npoints = -1 , int nfiles = -1 ) const;
     
     // mode 0 - File instanses are padded with 6 zeros
     // mode 1 - File instanses are not padded at all 
     void     SetMode  ( int x ) {fMode = x; }

     TString GetFileInstanceString     () { return FileInstance ;}
//   static TChain * CheckAndBuildChain       ( CardReader * MasterCard, std::string mode );
     static TChain * CheckAndBuildChain   ( CardReader * MasterCard, std::string mode, std::string instance="-1"  , std::string rep="*" );

     void UseNew(bool kNew=true){kUseNew=kNew;};

  private: 

     int       Instance;
     int       nPoints;
     int       nFiles;
     int       fMode ;
     TString   FilePath;
     TString   FileInstance;
     bool kUseNew;


};

#endif


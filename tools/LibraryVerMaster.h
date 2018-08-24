#ifndef __LibraryVerMaster__
#define __LibraryVerMaster__

#include <iostream>
#include <string>
#include <map>

#include "TFile.h"
#include "TTree.h"
#include "TString.h"

class LibraryVerMaster
{
 public :

   LibraryVerMaster();

   static LibraryVerMaster * GetMaster();

   void   SaveVersionInfo( );
   void   SaveVersionInfo( TFile * );

   void   AddLibrary( LibraryVerMaster * ptr );

   const char * GetName(){ return name.c_str(); }

   int GetVersion    () {return Ver;}
   int GetSubVersion () {return Sub;}
   int GetRevision   () {return Rev;}
 protected: 

   static LibraryVerMaster * theMaster ;

   std::string name;
   int Ver;
   int Sub;
   int Rev;

   static std::map< std::string , LibraryVerMaster * > regLibs;
 
};

#endif

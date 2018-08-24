#ifndef __coreLibVer__
#define __coreLibVer__

#include <iostream>
#include <string>

#include "tools/LibraryVerMaster.h"

class coreLibVer : public LibraryVerMaster
{
 public :

   coreLibVer();

   static coreLibVer * GetLibrary();

 protected: 

   static coreLibVer * theMaster ;
};

#endif

#ifndef __toolsLibVer__
#define __toolsLibVer__

#include <iostream>
#include <string>

#include "tools/LibraryVerMaster.h"

class toolsLibVer : public LibraryVerMaster
{
 public :

   toolsLibVer();

   static toolsLibVer * GetLibrary();

 protected: 

   static toolsLibVer * theMaster ;
};

#endif

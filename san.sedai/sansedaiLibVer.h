#ifndef __sansedaiLibVer__
#define __sansedaiLibVer__

#include <iostream>
#include <string>

#include "tools/LibraryVerMaster.h"

class sansedaiLibVer : public LibraryVerMaster
{
 public :

   sansedaiLibVer();

   static sansedaiLibVer * GetLibrary();

 protected: 

   static sansedaiLibVer * theMaster ;
};

#endif

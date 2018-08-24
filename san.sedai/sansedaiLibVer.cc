#include "san.sedai/sansedaiLibVer.h"

sansedaiLibVer * sansedaiLibVer::theMaster = 0;

sansedaiLibVer::sansedaiLibVer()
{  
  name = "sansedaiLibVer";
  Ver  =  2 ;
  Sub  =  0 ;
  Rev  =  0 ;
}

sansedaiLibVer * sansedaiLibVer::GetLibrary()
{
  if ( sansedaiLibVer::theMaster == NULL ) 
  {
    sansedaiLibVer::theMaster = new sansedaiLibVer;
    return sansedaiLibVer::theMaster;
  }
  else 
       return sansedaiLibVer::theMaster;
}


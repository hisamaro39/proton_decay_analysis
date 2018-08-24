#include "core/coreLibVer.h"

coreLibVer * coreLibVer::theMaster = 0;

coreLibVer::coreLibVer()
{  
  name = "coreLibVer";
  Ver  =  1 ;
  Sub  =  3 ;
  Rev  =  1 ;
}

coreLibVer * coreLibVer::GetLibrary()
{
  if ( coreLibVer::theMaster == NULL ) 
  {
    coreLibVer::theMaster = new coreLibVer;
    return coreLibVer::theMaster;
  }
  else 
       return coreLibVer::theMaster;
}


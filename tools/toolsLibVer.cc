#include "tools/toolsLibVer.h"

toolsLibVer * toolsLibVer::theMaster = 0;

toolsLibVer::toolsLibVer()
{  
  name = "toolsLibVer";
  Ver  =  1 ;
  Sub  =  0 ;
  Rev  =  1 ;
}

toolsLibVer * toolsLibVer::GetLibrary()
{
  if ( toolsLibVer::theMaster == NULL ) 
  {
    toolsLibVer::theMaster = new toolsLibVer;
    return toolsLibVer::theMaster;
  }
  else 
       return toolsLibVer::theMaster;
}


**********************************************************************
*     -------------------
      SUBROUTINE NECLRCRS
*     -------------------
*     
*     ( purpose )
*       Clear common NECOH
*
*     ( input ) 
*       none
*
*     ( output ) 
*       COMMON /NECOH/
*
*     ( creation date and author )
*       2011.5.23 ; First version by P. de Perio
*
*     ( comment )
*      
**********************************************************************

      IMPLICIT NONE

      integer*4 i

#include "neutcrs.h"

      CRSENERGY = -1
      CRSX = -1
      CRSY = -1
      CRSZ = -1
      CRSPHI = -1
      
      TOTCRSNE = -1
      DO i = 1, 8
         DIFCRSNE(i)=-1
      enddo

      RETURN
      END

#ifndef _AlignFriend_
#define _AlignFriend_

#include <iostream>
#include <string>
#include "DataManager.h"
#include "Rtypes.h"


class AlignKeyParameters
{
  public:
   AlignKeyParameters( const char * s , int bs, int nb , int sp )
   {
      Var           = s;
      BitStart      = bs;
      nBits         = nb;
      StuffPosition = sp;
   }

   std::string Var; 
   int         BitStart;
   int         nBits;
   int         StuffPosition;

}; 



class AlignFriend
{
  public:
    AlignFriend( DataManager * );

    void ClearInputKeys( ) { InputKeys.clear() ; }
    void AddInputKey( const char * s , int bs, int nb , int sp );

    ULong64_t BuildAlignmentVar();

  protected:
    virtual void ConsistencyTest() {};
    
    DataManager * dm; 
    std::vector< AlignKeyParameters * > InputKeys;

    void     DisplayBits( ULong64_t target  , unsigned bytes , const char * txt);
    unsigned StuffBits( ULong64_t & target , unsigned object , int start, int nbits );
    unsigned ExtractBits( unsigned target , int start, int nbits );

};



#endif

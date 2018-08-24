#include "tools/AlignFriend.h"

AlignFriend::AlignFriend( DataManager * ldm )
{
 
  dm = ldm ;

  ConsistencyTest();

}


void AlignFriend::AddInputKey(  const char * s , int ns , int nb , int sp )
{
  InputKeys.push_back( new AlignKeyParameters(  s , ns, nb , sp ) );

  // this line is here to force our parameters to be 
  // filled in the use_list and hence readout on calls 
  // to dm->GetEntry(); 
  dm->AddToUseList( s );
}



ULong64_t AlignFriend::BuildAlignmentVar()
{


 if ( InputKeys.size() < 1 ) 
 {
    std::cout << "AlignFriend::BuildAlignmentVar Error - now InputKeys specified. " << std::endl;
    std::cout << "    An alignment ID cannot be formed, returning 0  " << std::endl; 
    return 0;
 }

 void     * vptr;
 unsigned * uptr;
 std::string var ;

 unsigned copy;
 unsigned x ;

 // defined to be 8 bytes
 ULong64_t AlignmentKey = 0 ; 

 for ( unsigned i = 0 ; i < InputKeys.size() ; i ++ )
 {

   var             = InputKeys[i]->Var;
   int BitStart    = InputKeys[i]->BitStart;
   int nStuff      = InputKeys[i]->nBits;
   int StuffStart  = InputKeys[i]->StuffPosition;

   vptr = dm->Get( var.c_str() );

   uptr = (unsigned*) ( vptr );
   copy = (*uptr);
    
   x = ExtractBits( copy , BitStart , nStuff );

// std::cout << "Var: " << var << " " << copy << " x: " << x << std::endl;
   StuffBits( AlignmentKey , x , StuffStart, nStuff ); 
   
 }

// DisplayBits( AlignmentKey   ,  8, "Key"    );

 return AlignmentKey;

}


void AlignFriend::DisplayBits( ULong64_t target , unsigned bytes , const char * txt )
{

   unsigned mask ;
   unsigned x ;
   unsigned size = bytes*8;

   char * rep = new char [ bytes ] ;   
   memcpy ( rep , &target , bytes );

   std::cout << " " << size << " " << txt << " "  << target << ": "  ;
   for( int i = bytes - 1 ; i >= 0 ; i-- )
   {
      for( int j = 7 ; j >=0 ; j-- )
         putchar( '0' + (( rep[i] >> j ) & 1 ) ) ;
      std::cout << " " ; 
   }
   std::cout << std::endl;  
   delete rep;


}


unsigned AlignFriend::StuffBits( ULong64_t & target , unsigned object , int start, int nbits )
{

   unsigned mask = 0 ;
   if( abs(nbits) < 32 ) 
      mask = ( 1 << abs(nbits) ) - 1 ;
   else 
      mask = ( 0xFFFFFFFF );

   ULong64_t local = 0 ;
   local  = ( object & mask   );
   local  = ( local  & mask   );
   local  = ( local << start  );
   target = ( target | local  );

}

unsigned AlignFriend::ExtractBits( unsigned target , int start, int nbits )
{

   unsigned x ;
   unsigned mask;


   if( abs(nbits) < 32 ) 
      mask = ( 1 << abs(nbits) ) - 1 ;
   else 
      mask = ( 0xFFFFFFFF );

   if ( nbits < 0 ) target = (~target) + 1 ;

   x  = ( target >> start );
   x  = (  x & mask ); 
   
  return x;
} 





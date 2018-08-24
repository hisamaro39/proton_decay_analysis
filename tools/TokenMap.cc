#include "tools/TokenMap.h"
#include <iostream>
#include <stdlib.h>


TokenList::TokenList()
{
   fTokens.clear();
   fNtokens = 0 ;
} 

TokenList::TokenList( std::vector< std::string > & v )
{
   fTokens.clear();
   fNtokens = v.size() ;
 
   fTokens = v ;
} 

void TokenList::AddToken( const char * c )
{
   std::string s(c);
   AddToken( s );
}

void TokenList::AddToken( std::string s )
{
   fTokens.push_back(s);
   fNtokens++;
}






void TokenMap::AddTokenList( const char * key, int ntokens, std::vector<std::string>& v )
{

  std::string sKey ( key );

  if( mTokenListMap.count( sKey ) )
  {
    std::cout << " TokenMap::AddTokenList Warning: Token Key " << sKey << " already exists " 
              << " and will be overwritten. " << std::endl; 

    delete mTokenListMap[ sKey ];
  } 

  mTokenListMap[ sKey ] = new TokenList( v ); 

}

void TokenMap::GetListOfKeys( std::vector<std::string>& v )
{
   std::map< std::string, TokenList * >::iterator _i;

   v.clear();
   for( _i = mTokenListMap.begin() ; _i != mTokenListMap.end() ; _i++ )
      v.push_back( _i->first );  

}

TokenList *  TokenMap::GetTokenList( const char * c )
{
   std::string sKey(c);
   if( ! mTokenListMap.count( sKey ) )
   {
      std::cout << " TokenMap::GetTokenList Warning: Token Key " << sKey << " not found. " 
                << " returning 0. " << std::endl; 
      return 0;
   }

   return mTokenListMap[ sKey ];
}

unsigned TokenMap::GetNtokens( const char * key )
{

   TokenList * list =  GetTokenList( key ); 

   if( list != 0 ) return list->GetNtokens();
   else               return 0;

}

std::string TokenMap::GetToken( const char * key , unsigned i )
{

   TokenList * list =  GetTokenList( key ); 

   if( list != 0 ) return list->GetToken(i);
   else               return "\0";

}

void TokenMap::GetToken( const char * key , unsigned entry, const char * o , void * ptr  )
{

   std::string opt(o);
   TokenList * list =  GetTokenList( key ); 

   if( list == 0 ) 
   {
     ptr = 0;
     return;
   }

   std::string Token = list->GetToken(entry);
   ProcessToken( ptr , Token, opt );

   return ;
}


void TokenMap::ProcessToken( void * ptr ,std::string Token , std::string opt )
{
   if( opt == "i" || opt == "I" ){ int    * i = (int*   ) ptr ;  *i =         atoi( Token.c_str() );  }
   if( opt == "f" || opt == "F" ){ float  * f = (float* ) ptr ;  *f = (float) atof( Token.c_str() );  }
   if( opt == "d" || opt == "D" ){ double * d = (double*) ptr ;  *d =         atof( Token.c_str() );  }

   if( opt == "s" || opt == "S" ){ std::string * s = (std::string*) ptr ;  *s = Token;      }

}










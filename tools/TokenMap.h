#ifndef _TokenMap_
#define _TokenMap_

#include<map>
#include<vector>
#include<string>


class TokenList
{
   public:
     TokenList( );
     TokenList( std::vector<std::string> & );

     void AddToken( const char * );
     void AddToken( std::string  );
     
     std::string GetToken( unsigned i ) { return fTokens[i]; }
     unsigned GetNtokens() { return fNtokens; } 

   private:
     std::vector< std::string > fTokens;
     unsigned fNtokens;

};


class TokenMap 
{
  public:
    TokenMap() { fKeyRoot = "" ; }  
    TokenMap( const char * key ) { fKeyRoot = key ; }

    void        AddTokenList( const char * , int , std::vector<std::string> & );
    void        GetListOfKeys( std::vector<std::string>&  );
    TokenList * GetTokenList( const char * );
  
    std::string GetKeyRoot() { return fKeyRoot; }
    std::string GetToken( const char * key , unsigned i );

    void        GetToken( const char * , unsigned , const char * , void * );
    void        ProcessToken( void *  ,std::string , std::string );
    unsigned    GetNtokens( const char * key );
 
  private:
    std::string fKeyRoot;
    std::map< std::string, TokenList * > mTokenListMap; 

};

#endif

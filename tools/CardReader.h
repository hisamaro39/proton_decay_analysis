#ifndef _CardReader_
#define _CardReader_

#include "tools/TokenMap.h"

#include <cstdlib>
#include <map>
#include <vector>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>

#include "TString.h"


class CardReader
{
   public:
     
      CardReader( const char *  , const char * , bool k = false );
      CardReader( const char *  , bool k = false );
      CardReader( std::string , bool k = false );

      void ReadCard( std::string );

      void KeyFindFailure( std::string key , const char * type , void * ptr );
      bool GetKey( std::string, double & );
      bool GetKey( std::string, float  & );
      bool GetKey( std::string, int    & );
      bool GetKey( std::string, bool   & );
      bool GetKey( std::string, std::string & );
      bool GetKey( const char*, std::string & );
  
      bool FindKey( std::string );

      int ParseArray( std::string , const char * , double ** );
      unsigned GetNtokens( const char * );

      double * Get( std::string a ) { return bins[a]; } 
      double * GetBinsForType( int  );

      void BuildListOfStrings( const char * , std::map< std::string, std::string > & );
      void MakeStringTokens( TString * line, const char * D , std::vector< TString *> & data );

      TokenMap * BuildTokenMap( const char * skey, const char * delim  );
      void GetListOfKeys( std::vector< std::string >& v );

      void Erase( std::string & , const std::string & , std::string::size_type l0 = 0 );

      unsigned int Tokenize(const std::string & source,
                            const std::string & delimiters,
                            std::vector<std::string> & tokens);

   private:
      void BuildBinMap();  
      void BuildDimensionMap();

      
      void Replace( std::string & , const std::string & , 
                                    const std::string & ,
                                    std::string::size_type l0 = 0 );

      std::map< std::string , double >  Keys;
      std::map< std::string , double >::iterator iMap;

      std::map< std::string , std::string >  StringKeys;
      std::map< std::string , std::string >::iterator iSMap;

      std::map< std::string , double * > bins;
      std::map< int    , std::string   > bintype;
      std::map< std::string , double * >::iterator ibMap;
      
      std::vector< char > commentChars;
    
      bool kVerbosity;

};




#endif




#include "tools/CardReader.h"
#include "TObjString.h"
#include "TObjArray.h"
#include "TObject.h"
#include <iostream>

#include "tools/toolsLibVer.h"

CardReader::CardReader( const char * cname , const char * comments , bool k)
{
   kVerbosity = k ;

   const char * ptr = comments;
   for( unsigned i = 0 ; i < strlen(comments) ; i++ )
   {
      if( kVerbosity ) 
          std::cout << "Adding comment char: " << *ptr << std::endl;
      commentChars.push_back( (*ptr++) );  
      
   }

   std::string CardName( cname );
   ReadCard( CardName );

   toolsLibVer * tlvp = toolsLibVer::GetLibrary();
   LibraryVerMaster::GetMaster()->AddLibrary( tlvp ); 
}

CardReader::CardReader( const char * cname , bool k)
{

   kVerbosity = k ;
   std::string CardName( cname );
   commentChars.push_back('#');
   ReadCard( CardName );

   toolsLibVer * tlvp = toolsLibVer::GetLibrary();
   LibraryVerMaster::GetMaster()->AddLibrary( tlvp );
}

CardReader::CardReader( std::string CardName , bool k)
{
   kVerbosity = k ;
   commentChars.push_back('#');
   ReadCard( CardName );

   toolsLibVer * tlvp = toolsLibVer::GetLibrary();
   LibraryVerMaster::GetMaster()->AddLibrary( tlvp ); 
}

void CardReader::ReadCard( std::string CardName )
{

   std::ifstream InFile;
   InFile.open( CardName.c_str(), std::ifstream::in );
	
   if( !InFile.is_open()  )
   {
      std::cout << "CardReader::ReadCard " << CardName << " is not valid. " << std::endl; 
   }
   std::cout << "CardReader::ReadCard reading input from " << CardName << std::endl;

   std::string line;
   std::string Value;
   std::string Name;
   std::string::size_type loc;

   int i = 1;
   int LineCount = 0;
   int ntokens_space;
   int ntokens_comma;
   int ntokens_string;
   bool kComment;
   std::vector < std::string > tokens;      
   std::vector < std::string > tokensTmp;      


   // Reset the list of keys
   Keys.clear();
   while( !InFile.eof() )
   {  
      getline( InFile, line);

      LineCount ++; 
      if ( LineCount > 10000 ) 
      { 
         std::cerr << "CardReader::ReadCard Too many lines in card file, 10000" << std::endl;
         abort();
      }

      kComment = false;
      for( unsigned k = 0 ; k < commentChars.size() ; k++ ) 
      {
        loc = line.find( commentChars[k] , 0 );      
        if( loc == 0 || line.empty() )
        {
           if( kVerbosity )
              std::cout << " Omitting Line as Comment: " << line << std::endl;
           kComment = true ;
        }
      }
      if( kComment ) continue;

      // turn any tabs into spaces
      Replace( line, "\t", " " );

      // search for simple variable definitions 
      ntokens_comma  = Tokenize( line , "," , tokensTmp ); 
      ntokens_space  = Tokenize( line , " " , tokensTmp ); 

 
      tokens.clear();
      Name  = tokensTmp[0];
      Value = "";
      for( unsigned k = 1 ; k < ntokens_space ; k++ )
         if( tokensTmp[k].length() > 0 &&  tokensTmp[k] != " " )
            tokens.push_back( tokensTmp[k] );
       
      // second loop so we do not place a spurious comma at the 
      // of sequences with a tail of spaces      
      for( unsigned k = 0 ; k < tokens.size() ; k++ )
      {
         Value += tokens[k];
         if( ntokens_comma == 1 )
           if( k < tokens.size()-1)
             Value += ","; 
      }

      // remove extra spaces
      Erase( Name , " " );
      Erase( Value, " " );

      // we are just strings now so forget about quotes 
      Erase( Value, "\"" ); 
   
      // remove braces for array variables
      Erase( Value, "{" ); 
      Erase( Value, "}" );

      StringKeys[ Name ] = Value;
    
      // here be kludge, fix it .
      if( tokens.size() == 1 )
         Keys[ Name ] = atof( Value.c_str() );
         

   }// end of loop over lines in card file

   if( kVerbosity )
      for( iSMap = StringKeys.begin() ; iSMap != StringKeys.end() ; iSMap++ )
          std::cout << iSMap->first << " [" << iSMap->second <<"]"<< std::endl;

}

void CardReader::KeyFindFailure( std::string key , const char * type , void * ptr )
{
  std::cout << "CardReader::GetKey \"" << key << "\" for type \"" << type << "\" was not found."
            << " variable at address " << ptr << " is unchanged." << std::endl; 
}


bool CardReader::FindKey( std::string SearchKey )
{
   return ( Keys.count( SearchKey ) > 0 ? true : false );
}


bool CardReader::GetKey( std::string SearchKey, double & V )
{
   if ( FindKey(SearchKey) )
   {
       V = Keys[ SearchKey ];
       return  true;
   }

   if( kVerbosity ) KeyFindFailure( SearchKey , "double" ,  &V );
   return false;
}


bool CardReader::GetKey( std::string SearchKey, float  & V )
{
   if ( FindKey(SearchKey) )
   {
       V = (float) Keys[ SearchKey ];
       return  true;
   }

   if( kVerbosity ) KeyFindFailure( SearchKey , "float" ,  &V );
   return false;
}
bool CardReader::GetKey( std::string SearchKey, int    & V )
{
   if ( FindKey(SearchKey) )
   {
       V = (int) Keys[ SearchKey ];
       return  true;
   }

   if( kVerbosity ) KeyFindFailure( SearchKey , "int" ,  &V );
   return false;
}
bool CardReader::GetKey( std::string SearchKey, bool   & V )
{
   if ( FindKey(SearchKey) )
   {
       V = Keys[ SearchKey ] > 0 ? true : false;
       return  true;
   }

   if( kVerbosity ) KeyFindFailure( SearchKey , "bool" ,  &V );
   return false;
}

bool CardReader::GetKey( std::string SearchKey, std::string & V )
{

   if ( StringKeys.count(SearchKey) > 0  )
   {
       V = StringKeys[ SearchKey ];
       return  true;
   }

   if( kVerbosity ) KeyFindFailure( SearchKey , "string" ,  &V );
   return false;
}
 
bool CardReader::GetKey( const char * s , std::string & V )  
{
   std::string SearchKey( s );
   return GetKey( SearchKey , V );
}

/// will automatically cut spaces and {} out std::string
int CardReader::ParseArray( std::string Value, const char *  D, double ** Out )
{

   double dValue;
   double  * tmp;

   std::vector< std::string > tokens;
   
   int ntokens = Tokenize( Value , D , tokens );   

   // add two spaces so that [0] stores number of bins
   ntokens+= 2;

   tmp = new double [ ntokens ];
   tmp[0] = (double)  ntokens  ;

   for( unsigned int i = 0 ; i < tokens.size(); i ++ )
     tmp[i+1] = atof( tokens[i].c_str() ); 
     

   *Out = &tmp[0];

   // the number of elements found
   return ntokens-2;

}


void CardReader::GetListOfKeys( std::vector< std::string >& v )
{
   v.clear();
   for( iSMap = StringKeys.begin() ; iSMap != StringKeys.end() ; iSMap++ )
      v.push_back( iSMap->first );

   return;
}

unsigned CardReader::GetNtokens( const char * c )
{
   std::string key(c);

   if( StringKeys.count( key ) == 0 ) 
      return 0 ;   

   std::vector<std::string> v;
         
   unsigned nTokens = Tokenize( StringKeys[key] , ",", v ); 

   return nTokens; 
}



void CardReader::BuildListOfStrings( const char * skey , std::map< std::string, std::string > & v )
{

   std::string::size_type loc=0;  
   std::string searchkey( skey );
   std::string key;
   v.clear();

   for( iSMap  = StringKeys.begin() ; iSMap != StringKeys.end() ; iSMap++ )
   {
       key = iSMap->first;
       loc = key.find( searchkey.c_str(), 0, searchkey.length() );       

       // skip non matching std::strings
       if( loc == std::string::npos ) continue;

       v[ iSMap->first ] = iSMap->second ;
   }

  return;
}




void CardReader::BuildBinMap()
{

   std::string key;
   std::string binname;
   std::string::size_type loc=0;  

   double * edges;
   int nbins;

   int counter = 0;
   for( iSMap  = StringKeys.begin() ; iSMap != StringKeys.end() ; iSMap++ )
   {
       key = iSMap->first;
       loc = key.find("bin_", 0 , 4);       

       if( loc == std::string::npos ) continue;

       // expect the prefix to be at the beginning(duh)
       binname = key.substr( loc + 4 );

       //std::cout << "Key: " << key << " " << binname << std::endl;
       nbins = ParseArray( iSMap->second , "," , &edges ); 
       //std::cout <<  binname << " " << nbins<< " " << " " << edges[0] << " " << edges[1] << " " << edges[nbins] << std::endl;

       bins   [ binname       ] = &edges[0];        
       bintype[ int(edges[1]) ] = binname;        

       counter++;
   }


}


void CardReader::BuildDimensionMap()
{


   std::string key;
   std::string binname;
   std::string::size_type loc;  

   std::vector < std::string > tokens;
   std::vector < std::string > sub_tokens;
   int ntokens;
   int nsubtokens;

   double * edges;
   int nbins;

   int counter = 0;
   for( iSMap  = StringKeys.begin() ; iSMap != StringKeys.end() ; iSMap++ )
   {
       key = iSMap->first;
       loc = key.find("dimension_");       

       if( loc == std::string::npos ) continue;

       // expect the prefix to be at the beginning(duh)
       binname = key.substr( loc + 10 );

       // split into bin_type_1 (cut_type_1): bin_type_2 (cut_type_2)
       ntokens = Tokenize( iSMap->second , ":" , tokens );

       //simple tokenize
       for( int i = 0 ; i < ntokens ; i++ )
       {
          nsubtokens = Tokenize( tokens[i],"()", sub_tokens );    
          

       }

       counter++;
   }


}

void CardReader::Erase( std::string & source, const std::string& kill, std::string::size_type l0 )
{

   std::string::size_type loc0 = source.find( kill, l0 );
  
   while( loc0 != std::string::npos )
   {
      source.erase( loc0, 1);
      loc0 = source.find(kill, loc0);
   }


}



unsigned int CardReader::Tokenize(const std::string& source,
                                        const std::string& delimiters,
                                        std::vector<std::string>& tokens)
{
    std::string::size_type prev_loc0 = 0;
    std::string::size_type loc0 = 0;
    std::string::size_type loc1 = 0;
    std::string sub;
    unsigned int ntokens = 0;

    tokens.clear();

    loc0 = source.find_first_of(delimiters, loc0);
    while (loc0 != std::string::npos)
    {

        sub = source.substr(prev_loc0, loc0 - prev_loc0);

        tokens.push_back( sub ) ;
        ntokens++;


        loc0++;
        prev_loc0 = loc0;
        loc0 = source.find_first_of(delimiters, loc0);
    }


    if (prev_loc0 < source.length())
    {
        tokens.push_back(source.substr(prev_loc0));
        ntokens++;
    }

    return ntokens;
}




double * CardReader::GetBinsForType( int a )
{

   return bins[ bintype[a] ];


}


void CardReader::MakeStringTokens( TString * line, const char * D , std::vector< TString *> & data )
{
    unsigned i;
    for( i = 0 ; i < data.size() ; i++ )
       if( data[i] )
         delete data[i]; 
    data.clear();


    TObjArray * tokens = line->Tokenize(D);
    unsigned size = tokens->GetEntries();

    TObjString * s;
    for( i = 0 ; i < size ; i++ )
    {
       s = (TObjString*) tokens->At(i);
       if( s->GetString().Length() > 0 )
          data.push_back( new TString( s->GetString().ReplaceAll(" ",1) ) );
    }

}


void CardReader::Replace( std::string & source, const std::string & kill , 
                                                const std::string & rep  , std::string::size_type l0 )
{

   std::string::size_type loc0 = source.find( kill, l0 );
  
   while( loc0 != std::string::npos )
   {
      source.replace( loc0, kill.length(), rep );
      loc0 = source.find(kill, loc0);
   }

}


TokenMap * CardReader::BuildTokenMap( const char * skey, const char * delim  )
{
   
   std::map< std::string, std::string > aMap;
   std::map< std::string, std::string >::iterator _i;
   std::vector < std::string > tokens;      
   int ntokens;

   TokenMap * tokenMap = new TokenMap( skey );

   BuildListOfStrings( skey, aMap ); 
  
   for( _i = aMap.begin() ; _i != aMap.end() ; _i++ )
   {
      tokens.clear();
      ntokens = Tokenize( _i->second , delim , tokens );   
     
      tokenMap->AddTokenList( _i->first.c_str() , ntokens, tokens ); 
   } 


  return tokenMap;
}





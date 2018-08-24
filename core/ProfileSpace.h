#ifndef _ProfileSpace_
#define _ProfileSpace_

#include <iostream>
#include <fstream>
#include <string>
#include <map>

// Class defines a point in the Oscillation Profiling Space
//  aka the parameters used in the MNS matrix


class ProfileSpace
{
     public:

          ProfileSpace( int a, int b, double *c )
          {  
            tag = a; nDims = b; 
            Dims = new double [ nDims ];
            for( int i = 0 ; i < nDims; i++ )
              Dims[i] = c[i];   
                              
          }


          double GetProfileDim( int a ){  if( a < nDims ) return Dims[a]; else return 0; }
                      

          void Register( const char * name , int l , bool kVerbose = false )
          {
             if( l >= nDims )
             {
                std::cout << "ProfileSpace::Register error " 
                          << l << " > nDims [" << nDims << "]" << std::endl;
                return;
             }

             std::string entry( name );
             registry[ entry ] = l;
             if( kVerbose )
             std::cout << "ProfileSpace::Register : " << entry << " " << l << " defined " << std::endl;
          }

          double Get( const char * name )
          {   
             std::string entry( name );
             if( ! registry.count( entry ) ) 
             {
               std::cout << "ProfileSpace::Get warning " 
               << "registry [" << name << "] does not exist..." << std::endl;
               return 0.;
             }
             else 
               return Dims[ registry[entry] ];
          }


          void Print()
          {
            std::map<std::string, int>::iterator i;
            std::cout << "ProfileSpace::   tag " << tag  ;
            for( i = registry.begin() ; i != registry.end() ; i++ )
               std::cout << " " << i->first << " : " << Dims[ i->second ] ; 
               std::cout <<std::endl;
          }


          int GetTag()   { return tag; };

          ~ProfileSpace()
          {
            delete [] Dims;
          }



     private:
          int tag;
          int nDims;
          double *Dims;
          static std::map< std::string, int > registry;

};


#endif


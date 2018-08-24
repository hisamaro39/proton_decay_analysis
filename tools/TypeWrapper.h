#ifndef _TypeWrapper_
#define _TypeWrapper_

#include <iostream>
#include <cstdlib>
#include <string>

//
//  Wrapper in the sense that the T * buffer, buf,
//  points to something that already exits. That is,
//  we aren't creating any new storage for T objects
//
//  This class is intended to be used to provide a 
//  consistent interface between a variety of data types
//  so that, for instance, the same code can be used to 
//  handle data stored in either ROOT or PAW NTUPLE 
//  formats if the data are first loaded into 
//  TypeWrappers.
//

template <typename T>
class TypeWrapper
{
  public:

    TypeWrapper()
    {
       copy  = false;
       ndims = 0;
       bytes = 0;
    }

    TypeWrapper( const TypeWrapper<T> &b )
    {
        copy  = true;
        bytes =  b.bytes;
        buf   =  b.buf;
        ndims =  b.ndims;
        dims  = &b.dims[0];     
        Zero  =  b.Zero;
    }


    TypeWrapper( T* a , int d[], int n, std::string ln = "" )
    {
       copy  = false;

       //name  = ln  ;

       Zero  = new T;
      *Zero  = 0;

       bytes = 1;
       ndims = (n <= 0 ? 1 : n );
       dims  = new int [ndims];
 
       buf = a;
       //std::cout << "  TypeWrapper::TypeWrapper ndims " << ndims << std::endl;
       for( unsigned i = 0 ; i < ndims ; i++ )    
       {
          //std::cout << "  TypeWrapper::TypeWrapper " << d[i] << std::endl;
          dims[i] = d[i];
          bytes  *= d[i];  
       }
       bytes *= sizeof( T );

       // Print();
    }

   ~TypeWrapper()
   {
      if( !copy && ndims > 0 )
      {
        delete [] dims;
        delete Zero;
      }

   }

//  CINT chokes on this ?
//  T operator*(){ return buf[0]; }
    T operator()(){ return buf[0]; }
    T operator()(int c) 
      { 

         //std::cout << " TypeWrapper::operator() " << dims[0] << " " << buf[c] << std::endl;
         if( c >= dims[0] || ndims != 1) 
         {
            std::cout << "Warning TypeWrapper::operator( int ) "
                      << "requested dimension " << c << ", exceeds array length, [" << dims[0] <<"] "
                      << "or number of dimensions is incorrect, " << ndims  
                      << " (returning zero." << std::endl; 
            this->Print();
            abort();
            return *Zero;
         }

         return buf[c]; 
      }


    T operator()(int r, int c) 
    {
       if( ndims != 2 || r >= dims[0] || c >= dims[1] )
       {
          std::cout << "Warning TypeWrapper::operator( int, int ) "
             << "requested dimension(s) exceed allocated space. " 
             << "or number of dimensions is incorrect, " << ndims  
             << "returning zero." << std::endl; 
          this->Print();
          return *Zero;
       }

       return buf[ dims[1]*r + c ]; 
    }

    T operator()(int i, int j, int k) 
    {
       if( ndims != 3 || i >= dims[0] || j >= dims[1] || k >= dims[2] )
       {
          std::cout << "Warning TypeWrapper::operator( int, int, int ) "
             << "requested dimension(s) exceed allocated space. " 
             << "or number of dimensions is incorrect, " << ndims  
             << "returning zero." << std::endl; 
          this->Print();
          return *Zero;
       }

       return buf[ dims[2]*dims[1]*i + dims[2]*j + k ]; 
    }


    // four dimensions
    T operator()(int i, int j, int k , int l) 
    {
       if( ndims != 4 || i >= dims[0] || j >= dims[1] || k >= dims[2] || l >= dims[3]  )
       {
          std::cout << "Warning TypeWrapper::operator( int, int, int ) "
             << "requested dimension(s) exceed allocated space. " 
             << "or number of dimensions is incorrect, " << ndims  
             << "returning zero." << std::endl; 
          this->Print();
          return *Zero;
       }

       return buf[ dims[3]*dims[2]*dims[1]*i + dims[3]*dims[2]*j + dims[3]*k + l ]; 
    }

    // five dimensions
    T operator()(int i, int j, int k , int l, int m) 
    {
       if( ndims != 5 || i >= dims[0] || j >= dims[1] || k >= dims[2] 
                      || l >= dims[3] || m >= dims[4] )
       {
          std::cout << "Warning TypeWrapper::operator( int, int, int ) "
             << "requested dimension(s) exceed allocated space. " 
             << "or number of dimensions is incorrect, " << ndims  
             << "returning zero." << std::endl; 
          this->Print();
          return *Zero;
       }

       return buf[ dims[4]*dims[3]*dims[2]*dims[1]*i + 
                   dims[4]*dims[3]*dims[2]        *j + 
                   dims[4]*dims[3]                *k + 
                   dims[4]                        *l +  m  ]; 
    }

    //
    // slightly deeper than shallow copy...
    //  
    TypeWrapper<T>& operator=( const TypeWrapper<T> &b )
    {
        copy  = true;
//      name  =  b.name;
        buf   =  b.buf;
        bytes =  b.bytes;
        ndims =  b.ndims;
        dims  = &b.dims[0];     
        Zero  =  b.Zero;

        return *this;
    }


    void Print()
    {
    // std::cout << "Start Of TypeWrapper " << name << std::endl;
       std::cout << "Start Of TypeWrapper " << std::endl;
       std::cout << "Type Size: " << sizeof( T ) << std::endl;
       std::cout << "copy:  " << copy  << " [ 0: not a copy, 1: is a copy ] " << std::endl;
       std::cout << "bytes: " << bytes << std::endl;
       std::cout << "ndims: " << ndims << std::endl;
       std::cout << "  dimensions: " << std::endl;
       for( int i = 0 ; i < ndims ; i++ )
          std::cout << "     dims[" << i << "] : " << dims[i] << std::endl;

       int max =1;
       for( int i = 0 ; i < ndims ; i++ )
          max *= dims[0];

       std::cout << std::endl;        
       std::cout << "buffer start: " << buf << std::endl;
       std::cout << "buffer contents: " << std::endl;
       std::cout << "   ";
       for( int i = 0 ; i < max ; i++ )
       {
          std::cout << buf[i] << " "; 
          if( i + 1 % 10 == 0 )
            std::cout << std::endl << "   ";
       }
       std::cout << std::endl;
       std::cout << std::endl;
       std::cout << "End Of TypeWrapper " << std::endl;
       std::cout << std::endl;

    }


  private:

    int    bytes;
    int    ndims;
    int *  dims;

    // assumes <T> is laid out lineraly in memory
    T *    buf;
 
    // were we made from copy constructor
    bool   copy; 

    T *   Zero;
    
//  std::string name;
};






#endif

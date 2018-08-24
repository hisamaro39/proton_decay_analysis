
#ifndef _Fortran_
#define _Fortran_


namespace Fortran
{

// deprecated
//extern "C" double upmuflx_( int & , int & , int & );

// neutrino path lenghts including production height from honda flux
extern "C" void  nebaseline_h3d_  ( float [], int& , float& , float&, float& ); 

// neutrino production heights
extern "C" void  neprodhgt_h3d_   ( float [], int& , float& , float&  ); 

// interface to fortran systematic error computation
// deprecated
//extern "C" double func_cf_( int&, int& , int&, int& );

// corresponding 1-sigma values 
// deprecated
// extern "C" double fgetsys_( int& );

}


#endif

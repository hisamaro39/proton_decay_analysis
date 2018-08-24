#ifndef _GenericOscPoints_
#define _GenericOscPoints_

#include <string>
#include <iostream>

class GenericOscPoints {

 public:
   GenericOscPoints( double start, double stop, int islog , int n , const char * name,  int wrappable = 0 )
   {
      Start      = start;  
      Stop       = stop; 
      kIsLog     = islog; 
      nSamples   = n; 
      Name       = name;
      PointList  = 0;
      Edges      = NULL;
      kWrappable = wrappable;
      CreatePointList();
   }


   GenericOscPoints( int npoints , const char * name, int wrappable = 0 )
   {
      kIsLog     = 2;
      nSamples   = npoints;
      Name       = name ;
      Start      = 0.0;
      Stop       = 0.0;
      Edges      = NULL;
      kWrappable = wrappable;
      PointList  = new double [ nSamples ];
   }

   ~GenericOscPoints( )
   {
     if ( PointList ) delete PointList ;
     if ( Edges )     delete Edges;
   }


   void SetPoint( unsigned i , double val )
   {
     if ( !PointList )
     {
        std::cout << "GenericOscPoints::SetPoint Error. PointList has not been created for parameter type " << Name<< std::endl;
        return;
     } 
     if( i < nSamples )
       PointList[i] = val;
     else 
       std::cout << "GenericOscPoints::SetPoint Error. Specified point " << i << " is larger than expected PointList size, "<< nSamples << std::endl;
   }


   double GetPoint( unsigned i )
   {
     if ( !PointList )
     {
        std::cout << "GenericOscPoints::GetPoint Error. PointList has not been created for parameter type " << Name << std::endl;
        return 0.0;
     } 
     if( i < nSamples )
       return PointList[i];
     else 
     {
       std::cout << "GenericOscPoints::SetPoint Error. Specified point " << i << " is larger than expected PointList size, "<< nSamples << std::endl;
       return 0.0;
     }
   }


   std::string Name;
   double Start;
   double Stop;
   unsigned int kIsLog;
   unsigned int nSamples;
   unsigned int kWrappable;

   double GetLowEdge( unsigned i ) 
   {
     if ( !Edges ) this->SetEdges();

     if( i <= nSamples + 1 )
       return Edges[i] ;
     return 0.0;
   }

   void SetEdges() 
   {
     if ( Edges ) delete Edges;
   
     Edges = new double [ nSamples + 2 ];

     double half;
     double StepSize;
     for ( unsigned i = 0 ; i < nSamples - 1 ; i++ )
     {   
        StepSize = ( PointList[i+1] - PointList[i] ); 
        if ( kIsLog == 1 )
          StepSize = log10( PointList[i+1] / PointList[i] ); 
        StepSize *= 0.5;

        // adjust bin edges to center around chosen osc points
        half = ( kIsLog == 1 ? PointList[i] * pow (10., StepSize ) 
                             : PointList[i] + StepSize ); 
        Edges[i+1] = half;
        if( i == 0 )
        {
          half = ( kIsLog == 1 ? PointList[i] / pow (10., StepSize ) 
                               : PointList[i] - StepSize ); 
          Edges[0] = half;
        }
     }
     // Use previous stepsize for upper edges
     unsigned index =  nSamples - 1;
     half = ( kIsLog == 1 ? PointList[ index ] * pow (10., StepSize ) 
             : PointList[index] + StepSize ); 
     Edges[ index + 1]  = half ;

     half = ( kIsLog == 1 ? PointList[ index ] * pow (10., 2. * StepSize ) 
             : PointList[index] + 2.0 * StepSize ); 
     Edges[ index + 2 ]  = half ;
   }
   
 protected:

   double * PointList;
   double * Edges;

   void CreatePointList()
   {
     PointList = new double [ nSamples ];

     double lnPoints  = (double) nSamples;

     double localValue = 0.0;
     if ( nSamples > 1 )
       for ( unsigned i = 0 ; i < nSamples ; i++ )
       {
         double j = double (i) ;

         if ( kIsLog == 1 )
          localValue = Start * pow( Stop/Start , j / (lnPoints-1.0) );

         if ( kIsLog == 0 )
          localValue = Start + ( Stop - Start ) * j / (lnPoints-1.0);
         SetPoint( i , localValue );
       }
     else 
       SetPoint( 0 , Start );
 
   }

};


#endif

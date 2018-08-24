

#ifndef _ThreeFlavorOscillator_
#define _ThreeFlavorOscillator_

#include "core/Oscillator.h"
#include "skcore/SKEventParser.h"
#include "Prob3++/BargerPropagator.h"
//#include "LineIntegral.h"
#include "tools/CardReader.h"


class ProfileSpace;
class EventParser;

class ThreeFlavorOscillator : public Oscillator
{
	public:

		//ThreeFlavorOscillator(bool );
		//ThreeFlavorOscillator( CardReader * );
		
		//double operator()( EventParser* , ProfileSpace& , int );

		//double Oscillate      ( int,  EventParser* , int  , ProfileSpace&  );
		double Oscillate      ( TTree*  );
                //double PathAve3D      ( int , EventParser* , int  , ProfileSpace&  );
                //double Official2D     ( int , EventParser* , int  , ProfileSpace&  );
                //double NoAve2D        ( int , EventParser* , int  , ProfileSpace&  );
                //double CPT2D          ( int , EventParser* , int  , ProfileSpace&  );
                //double CPT2DAngleFix  ( int , EventParser* , int  , ProfileSpace&  );
		//double ENWmodel       ( int , EventParser* , int  , ProfileSpace&  );
                //double Neighbor3D     ( int , EventParser* , int  , ProfileSpace&  );

 	private:
                //BargerPropagator *bNu;

		//int nPaths;

		//int OscFlag;
                //bool kInverted;

		//float lContribution;
};
#endif

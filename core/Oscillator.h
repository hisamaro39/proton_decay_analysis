

#ifndef _Oscillator_
#define _Oscillator_

#include "core/EventParser.h"
#include "core/profiler.h"
#include "core/ProfileSpace.h"
#include <iostream>
#include <math.h>

class ProfileSpace;
class EventParser;

class Oscillator
{
	public:
		Oscillator(){};

		virtual double Oscillate( int , EventParser* , int, ProfileSpace& )=0;

	protected:


};
#endif

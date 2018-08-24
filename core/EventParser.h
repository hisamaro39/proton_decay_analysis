#ifndef _EventParser_
#define _EventParser_
/*  
   The purpose of the EventParsing class is to extract the variables relevant to an analysis
   from a data file and collect them, one event at a time. In many instances of the code, 
   EventParsers are used as one might use the event itself. For SK there are internal variables,
   disctinct from those in the event ntuples, which are used in the making of the Systabe and 
   are gathered, and accessed in this class.

   The User should derive a class and override the virtual functions, based on the details of their
   analysis. For T2K, as an example, the event structure is much simpler than that of SK, so the 
   parse routine should be overriden to accommodate these differences.

*/

#include "tools/DataManager.h"
#include "core/Event.h"
#include <math.h>
#include <iostream>
#include <vector>

class EventParser
{
  public:

    EventParser();
    virtual ~EventParser() {;};
    
    int GetPDG() 	    { return PDG; };		/* Primary particle ID at Event vertex in PDG	*/
    int GetMode()	    { return Mode; };		/* Interaction code at the vertex		*/

    virtual void Print()    { std::cout << " Event Type: " << EventType << std::endl; }

    double GetMCWeight()    { return MCweight; };
    double GetDataWeight()  { return DTweight; };
  
    int GetParserType()     { return ParserType ; }


    bool MCKaNa()	    { return MonteCarlo;};
    bool Skip()	            { return SkipFlag;};
    /* Should this Event be skipped for a given analysis- good for selecting events
       out of larger samples									*/

    virtual int GetType()   { return EventType; }; 

    virtual double GetOscProb() {return 0;};
    /* In case of neutrino analysis the Oscillation probability should be redefined			*/

    virtual void   Parse() = 0 ;
    virtual void   GetEventBinValue ( int ) {};
    virtual void   GetEventBinValues( std::vector<double> & ) {};
    virtual void   LoadDataManager( DataManager * ) = 0 ;

    void  SetTypeOffset ( int x ) { TypeOffset = x ; }


    void   RegisterVar( std::string , void* );
    void * GetVar     ( std::string );

  protected:

    int   lun;		/* Fortran File numbering system  */

    int    EventType;	/* EventType 		*/
    int    PDG;        /* PDG code 								*/
    int    Mode; 		/* interaction type 							*/
    
    bool   SkipFlag;        /* Skip this event? 							*/
    double MCweight;
    double DTweight;

    bool   MonteCarlo;
    int    ParserType ;
   
    int    TypeOffset ;

    std::map<std::string , void * > VarMap; 

};

#endif

#include "core/EventParser.h"

EventParser::EventParser()
{
   EventType = 0 ; 
   SkipFlag  = false ; 
}

void EventParser::RegisterVar( std::string s , void * p )
{
   VarMap[s] = p;
}

void * EventParser::GetVar( std::string s )
{
   int size;
   if( VarMap.count( s ) )
      return VarMap[s];
   else
      std::cout << " Erorr in EventParser::GetVar , - variable -" << s << "- was not found. " << std::endl;

   return 0;
}


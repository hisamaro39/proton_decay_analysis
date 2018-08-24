#ifndef __llHistogram__
#define __llHistogram__


#include "TFile.h"
#include "TH1D.h"

#include <string>
#include <sstream>

class llHistogram
{

 public:

   llHistogram( int _id, int _nbins, float lowEdge , float hiEdge, const char * title = "" , const char * _var = "")
   {
   
      ID = _id;
      var = _var;

      std::stringstream ss;
      ss << "h" << ID ;

      h = new TH1D( ss.str().c_str() , title , _nbins, lowEdge, hiEdge );
      std::cout << "llHistogram: " << ss.str() << " nBins:" << _nbins << " Title: " << title << " created. " << std::endl;
   }
  
   int ID;
   std::string var;
   TH1D * h ;

};


#endif


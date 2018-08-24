#ifndef _global_
#define _global_ 

#include <string>
namespace global
{


	enum BinDimensions{ _PDim, _ZDim };
        
        enum FileTypes{ fSK1 , fSK2, fSK3, fSK4 };

        enum gSKx{ SK2 = 19 , SK3 = 38, SK4 = 57};
        enum ntype{ nEventTypes = 19 };

        enum Parsers{ skP };

	enum BinTypes{    
                          SubGeV_elike_0dcy,
                          SubGeV_elike_1dcy,
                          SubGeV_SingleRing_pi0like,
                          SubGeV_mulike_0dcy,
                          SubGeV_mulike_1dcy,
                          SubGeV_mulike_2dcy,
                          SubGeV_pi0like,
		          MultiGeV_elike_nue,
		          MultiGeV_elike_nuebar,
		          MultiGeV_mulike,
		          MultiRing_elike_nue,
 		          MultiRing_elike_nuebar,
		          MultiRing_mulike,
        	          MultiRingOther_1,
		          PCStop,
			  PCThru,
		          UpStop_mu,
                          UpThruNonShower_mu,
                          UpThruShower_mu,

                          SK2_SubGeV_elike_0dcy,
                          SK2_SubGeV_elike_1dcy,
                          SK2_SubGeV_SingleRing_pi0like,
                          SK2_SubGeV_mulike_0dcy,
                          SK2_SubGeV_mulike_1dcy,
                          SK2_SubGeV_mulike_2dcy,
                          SK2_SubGeV_pi0like,
		          SK2_MultiGeV_elike_nue,
		          SK2_MultiGeV_elike_nuebar,
		          SK2_MultiGeV_mulike,
		          SK2_MultiRing_elike_nue,
		          SK2_MultiRing_elike_nuebar,
		          SK2_MultiRing_mulike,
        	          SK2_MultiRingOther_1,
		          SK2_PCStop,
			  SK2_PCThru,
		          SK2_UpStop_mu,
                          SK2_UpThruNonShower_mu,
                          SK2_UpThruShower_mu,

                          SK3_SubGeV_elike_0dcy,
                          SK3_SubGeV_elike_1dcy,
                          SK3_SubGeV_SingleRing_pi0like,
                          SK3_SubGeV_mulike_0dcy,
                          SK3_SubGeV_mulike_1dcy,
                          SK3_SubGeV_mulike_2dcy,
                          SK3_SubGeV_pi0like,
                          SK3_MultiGeV_elike_nue,
                          SK3_MultiGeV_elike_nuebar,
                          SK3_MultiGeV_mulike,
                          SK3_MultiRing_elike_nue,
                          SK3_MultiRing_elike_nuebar,
                          SK3_MultiRing_mulike,
        	          SK3_MultiRingOther_1,
                          SK3_PCStop,
                          SK3_PCThru,
                          SK3_UpStop_mu,
                          SK3_UpThruNonShower_mu,
                          SK3_UpThruShower_mu,

                          SK4_SubGeV_elike_0dcy,
                          SK4_SubGeV_elike_1dcy,
                          SK4_SubGeV_SingleRing_pi0like,
                          SK4_SubGeV_mulike_0dcy,
                          SK4_SubGeV_mulike_1dcy,
                          SK4_SubGeV_mulike_2dcy,
                          SK4_SubGeV_pi0like,
                          SK4_MultiGeV_elike_nue,
                          SK4_MultiGeV_elike_nuebar,
                          SK4_MultiGeV_mulike,
                          SK4_MultiRing_elike_nue,
                          SK4_MultiRing_elike_nuebar,
                          SK4_MultiRing_mulike,
        	          SK4_MultiRingOther_1,
                          SK4_PCStop,
                          SK4_PCThru,
                          SK4_UpStop_mu,
                          SK4_UpThruNonShower_mu,
                          SK4_UpThruShower_mu,
			  EndOfBinTypes,
		      };


	extern std::string BinStrings [];

        extern double UpStop_BG [4][2][3]  ;
        extern double UpShower_BG[4]    ;
        extern double UpNonShower_BG[4] ;

}
#endif


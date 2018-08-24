#include "skcore/global.h"

#include <string>

// updated for 201606 SK-IV data
      double global::UpNonShower_BG[4] = { 16.6 , 16.7 , 4.1 , 13.3 };
      double global::UpShower_BG[4]    = { 13.3 , 20.9 , 8.0 , 16.8 };

                                           //     z1           z0
                                           //p0,  p1,  p2,   p0,   p1, p2 
      double global::UpStop_BG [4][2][3] = {0.4, 0.1, 0.0, 17.9,  6.9, 4.7,     //sk1   
                                            0.1, 0.1, 0.0,  6.9,  6.8, 1.7,     //sk2
                                            0.3, 0.0, 0.0, 10.1,  3.1, 2.8,     //sk3
                                            0.0, 0.0, 0.9, 34.0,  6.1, 3.7 };   //sk4


      std::string global::BinStrings[] = {
                         "SubGeV_elike_0dcy",
                         "SubGeV_elike_1dcy",
                         "SubGeV_SingleRing_pi0like",
                         "SubGeV_mulike_0dcy",
                         "SubGeV_mulike_1dcy",
                         "SubGeV_mulike_2dcy",
                         "SubGeV_pi0like",
		         "MultiGeV_elike_nue",
		         "MultiGeV_elike_nuebar",
		         "MultiGeV_mulike",
		         "MultiRing_elike_nue",
 		         "MultiRing_elike_nuebar",
		         "MultiRing_mulike",
        	         "MultiRingOther_1",
		         "PCStop",
			 "PCThru",
		         "UpStop_mu",
                         "UpThruNonShower_mu",
                         "UpThruShower_mu",

                         "SK2_SubGeV_elike_0dcy",
                         "SK2_SubGeV_elike_1dcy",
                         "SK2_SubGeV_SingleRing_pi0like",
                         "SK2_SubGeV_mulike_0dcy",
                         "SK2_SubGeV_mulike_1dcy",
                         "SK2_SubGeV_mulike_2dcy",
                         "SK2_SubGeV_pi0like",
		         "SK2_MultiGeV_elike_nue",
		         "SK2_MultiGeV_elike_nuebar",
		         "SK2_MultiGeV_mulike",
		         "SK2_MultiRing_elike_nue",
		         "SK2_MultiRing_elike_nuebar",
		         "SK2_MultiRing_mulike",
        	         "SK2_MultiRingOther_1",
		         "SK2_PCStop",
			 "SK2_PCThru",
		         "SK2_UpStop_mu",
                         "SK2_UpThruNonShower_mu",
                         "SK2_UpThruShower_mu",

                         "SK3_SubGeV_elike_0dcy",
                         "SK3_SubGeV_elike_1dcy",
                         "SK3_SubGeV_SingleRing_pi0like",
                         "SK3_SubGeV_mulike_0dcy",
                         "SK3_SubGeV_mulike_1dcy",
                         "SK3_SubGeV_mulike_2dcy",
                         "SK3_SubGeV_pi0like",
                         "SK3_MultiGeV_elike_nue",
                         "SK3_MultiGeV_elike_nuebar",
                         "SK3_MultiGeV_mulike",
                         "SK3_MultiRing_elike_nue",
                         "SK3_MultiRing_elike_nuebar",
                         "SK3_MultiRing_mulike",
        	         "SK3_MultiRingOther_1",
                         "SK3_PCStop",
                         "SK3_PCThru",
                         "SK3_UpStop_mu",
                         "SK3_UpThruNonShower_mu",
                         "SK3_UpThruShower_mu",

                         "SK4_SubGeV_elike_0dcy",
                         "SK4_SubGeV_elike_1dcy",
                         "SK4_SubGeV_SingleRing_pi0like",
                         "SK4_SubGeV_mulike_0dcy",
                         "SK4_SubGeV_mulike_1dcy",
                         "SK4_SubGeV_mulike_2dcy",
                         "SK4_SubGeV_pi0like",
                         "SK4_MultiGeV_elike_nue",
                         "SK4_MultiGeV_elike_nuebar",
                         "SK4_MultiGeV_mulike",
                         "SK4_MultiRing_elike_nue",
                         "SK4_MultiRing_elike_nuebar",
                         "SK4_MultiRing_mulike",
        	         "SK4_MultiRingOther_1",
                         "SK4_PCStop",
                         "SK4_PCThru",
                         "SK4_UpStop_mu",
                         "SK4_UpThruNonShower_m",
                         "SK4_UpThruShower_m",
			 "EndOfBinTypes"
		      };


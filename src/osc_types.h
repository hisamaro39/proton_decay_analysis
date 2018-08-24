#ifndef _osc_types_
#define _osc_types_


   enum osc_types {  
            // default to fortran indexing here 
            // so conversion to paw ntuple later 
            // is smooth
            SubGeV_elike_0dcy = 1,
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
            EndOfTypes,

//
//          MultiRingOther_1,
            MultiRingOther_2,
           };


    enum SKGenerations { SK1, SK2, SK3, SK4 };

#endif

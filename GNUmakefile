#
#  In case of building locally, set ATMPD_ROOT variable 
#      setenv ATMPD_ROOT   ... 
#  or directly set in this makefile 
# ATMPD_ROOT = /skofl
#
ifndef ATMPD_ROOT
  ATMPD_ROOT = ../../..
endif


include $(ATMPD_ROOT)/config.gmk

#osctools     = $(ATMPD_ROOT)/src/analysis/Osc3++/tools
#osccore      = $(ATMPD_ROOT)/src/analysis/Osc3++/core
#ll           = $(ATMPD_ROOT)/src/analysis/EventSelectionLikelihoods/


#LOCAL_INC       = -I. -I$(osctools) -I$(ll) -I./EventSelectionLikelihoods/ -I$(osccore)
LOCAL_INC       = -I.  

#LOCAL_LIBS      =  -lAccessProb -lThreeProb -losctools -losc3ppcore -lneutflux $(ll)/libeventselectionll.a
#LOCAL_LIBS      =  -lAccessProb -lThreeProb  -losc3ppcore -lneutflux ./ESL/libeventselectionll.a
LOCAL_LIBS      =  -lAccessProb -lThreeProb  -losc3ppcore ./neutflux/libneutflux_5.4.0.a ./ESL/libeventselectionll.a



#  Objects

OBJS = src/OscNtupleManager.o tools/LibraryVerMaster.o tools/toolsLibVer.o tools/TokenMap.o tools/CardReader.o tools/AlignFriend.o tools/DataManager.o core/FileRecord.o ESL/LikelihoodReturn.o ESL/LikelihoodHelper.o ESL/mgmreLikelihood.o ESL/Pi0Likelihood.o ESL/nuebarLikelihood.o ESL/numuLikelihood.o ESL/numu1RLikelihood.o ESL/pcstopLikelihood.o ESL/pcthruLikelihood.o ESL/LikelihoodManager.o tools/HistogramService.o tools/LibraryVerMaster.o san.sedai/sansedaiLibVer.o core/EventParser.o core/ProfileSpace.o core/coreLibVer.o core/profiler.o core/Oscillator.o skcore/global.o skcore/SKEventParser.o Prob3++/EarthDensity.o Prob3++/mosc.o Prob3++/mosc3.o Prob3++/BargerPropagator.o san.sedai/ThreeFlavorOscillator.o fortran/fort_func.o 



build_osc_ntuple: build_osc_ntuple.o $(OBJS) $(FOBJS)
	LD_RUN_PATH=$(LIBDIR) $(CXX) $(CXXFLAGS) -o build_osc_ntuple build_osc_ntuple.o $(OBJS) $(FOBJS) $(LDLIBS)


LIBNAME = OscNtuple

.PHONY:  lib$(LIBNAME).a $(LIBDIR)lib$(LIBNAME).a

lib$(LIBNAME).a : $(OBJS)
	$(RM) $@
	$(AR) $@ $(OBJS)
	$(RANLIB) $@

$(LIBDIR)lib$(LIBNAME).a : lib$(LIBNAME).a
	$(RM) $@
	$(INSTALL_LIB) $< $@



.cc.so: 
	$(CXX) $(CXXFLAGS) -c $*.cc


emptyrule:: lib

lib:: lib$(LIBNAME).a $(OBJS)

install.lib:: $(LIBDIR)lib$(LIBNAME).a


clean:
	rm -f $(OBJS) lib$(LIBNAME).a  build_osc_ntuple *.o 

.PHONY:  clean setup includes install.includes depend lib install.lib exec install.exec


emptyrule:: lib


setup::

includes:: $(INCFILES) 

depend::


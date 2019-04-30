#
#  In case of building locally, set ATMPD_ROOT variable 
#      setenv ATMPD_ROOT   ... 
#  or directly set in this makefile 
# ATMPD_ROOT = /home/atmpd/skrep/08a
#
ifndef ATMPD_ROOT
  ATMPD_ROOT = ../../..
endif

include $(ATMPD_ROOT)/config.gmk

#
#  Objects
#

NEUTFLXVER = 5.4.0

LIBBASE = neutflux
LIBNAME = neutflux_${NEUTFLXVER}

FLO = 	fnflufx.o fnhonfx03.o fnhonfx03k.o fnhonfx06.o neflux.o \
	fnflufx03.o fnhonfx_1d.o neflxgai.o nenauflx.o \
	fngaifx.o neflxhon.o neprodhgt_h3d.o \
	nernneut.o fnbarfx.o fnhdfxh_1d.o neflxhon_1d.o \
	nerdbes.o nernnudg.o fnhonfx11low.o fnhonfx11.o \
	fnbarfx03.o fnhdfxl.o negaiflx.o nerdflu.o nernnudh.o \
	fnbesfx.o fnhdfxl_1d.o negvcd.o nerdgai.o nernnudh_1d.o \
	fnbsfxl.o fnhdfxm.o nebaseline_h3d.o nehonflx.o nerdhon.o \
	fnbsfxm.o fnhdfxm_1d.o necardvc.o neleeflx.o nerdhon_1d.o \
	fnfdfx.o fnhonfx.o neevrate.o nemake.o \
	nerdlee.o nebaseline_h3d2.o nebaseline_h3d2l.o \
	fnfidmas.o fnhonfx01.o nefint.o nemkevnt.o nerdnau.o rnneue.o \
	necardap.o nevectap.o \
	rnpos_wtrtnk.o rnpos_wtr.o
#	nemkfsibk.F nemkmodelbk.F nemknebk.F nerdfsibk.F 
UPM =   vector3.o chk_mu_upmu.o chk_nu_upmu.o dirtrans.o hrange.o dxde.o \
	rnpos_rck.o

ifdef SITE_HAS_NEUT
PROGRAMS = neut
else
PROGRAMS =
endif

ifdef UPMU_RCK
FORTRANDEFINES += -DUPMU_RCK -DUPMU
endif
ifdef UseG77
FORTRANDEFINES += -ffortran-bounds-check -g -O0 
#FORTRANDEFINES += -ffortran-bounds-check
else
FORTRANDEFINES += -fbounds-check -g -O0 
endif

ifdef UPMU
OBJS = $(FLO) $(UPM)
else
OBJS = $(FLO)
endif



all:: lib exec
#	flux


#
#  Rules for building library 
#

.PHONY:  lib$(LIBNAME).a $(LIBDIR)lib$(LIBNAME).a

lib$(LIBNAME).a : $(OBJS)
	$(RM) $@
	$(AR) $@ $(OBJS) 
	$(RANLIB) $@

$(LIBDIR)lib$(LIBNAME).a : lib$(LIBNAME).a
	$(RM) $@
	$(INSTALL_LIB) $< $@
	( cd $(LIBDIR) ; $(RM) lib$(LIBBASE).a ; $(LN) lib$(LIBNAME).a lib$(LIBBASE).a )

#
#  Targets 
#

.PHONY:  clean setup includes install.includes depend lib install.lib exec install.exec


emptyrule:: lib


clean::
	$(RM) *.o *~ *.a *.so core* neut neutntpl

setup::

depend::

lib:: lib$(LIBNAME).a

install.lib:: $(LIBDIR)lib$(LIBNAME).a

exec:: $(PROGRAMS)

install.exec:: 

#
# local targets
#

INCFILES = gaisser.h hondaflx.h bessflx.h flukaflx.h \
           necardap.h 

CONSTFILES = bartol03.dat \
            bess71.dat bess74.dat \
            fit_fluka01.dat fit_fluka02.dat fit_fluka02_mx.dat \
            fluka00.dat fluka01.dat fluka03.dat \
            gaisser96.dat \
            hkkm03mt.dat hkkm06mt.dat hkkm2.dat hkkm11mt.dat hkkm11low.dat\
            honda96high.dat honda96low.dat honda96mid.dat honda97mid.dat \
            kam_phi.smin kam_phi.smax

POBJS = trapfpe.o

CLIB = `cernlib pdflib804 jetset74 photos202 mathlib packlib kernlib`

MACHINEDIR=${MACHINE}

neut : neut.o $(POBJS) 
	if test -f ${NEUT_ROOT}/lib/libneutcore_${NEUTCOREVER}.a ; then LD_RUN_PATH=$(LIBDIR):$(SKOFL_LIBDIR) $(CXX) $(CXXFLAGS) -o neut $@.o $(POBJS) $(OBJS) ${NEUTLIBS} $(LDLIBS) ; fi

neutntpl : neutntpl.o $(POBJS)
	LD_RUN_PATH=$(LIBDIR):$(SKOFL_LIBDIR) $(CXX) $(CXXFLAGS) -o neutntpl $@.o $(POBJS) $(OBJS) ${NEUTLIBS} $(LDLIBS) 

includes:: $(INCFILES) 

install.includes:: $(INCFILES)
	$(INSTALL_INC) -t $(FORTRAN_INCDIR) $(INCFILES) 

install.lib:: $(CNSTFILES)
	$(INSTALL_DAT) -t $(CONSTDIR) $(CONSTFILES) 

clean::
	list='$(INCFILES)'; for file in $$list; do \
	( $(RM) $(FORTRAN_INCDIR)/"$$file" ); \
	done
	list='$(CONSTFILES)'; for file in $$list; do \
	( $(RM) $(CONSTDIR)/"$$file" ); \
	done


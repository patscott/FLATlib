# Makefile for CUBPACK within flatlib
#
# Pat Scott Feb 2009
# pat@fysik.su.se

LOCAL_ROOT = $(PWD)
LOCAL_BUILD = ${LOCAL_ROOT}/build
LOCAL_SRC = ${LOCAL_ROOT}/src

CORE_OBJ_BARE = buckley.o internal_types.o ds_routines.o divide.o \
	rule_tn.o rule_t3.o rule_t2.o rule_c2.o rule_c3.o rule_cn.o  \
	rule_1.o rule_general.o region_processor.o volume.o \
	check.o global_all.o error_handling.o cui.o
CORE_OBJ = $(CORE_OBJ_BARE:%.o=$(LOCAL_BUILD)/Core/%.o)

DRIVERS_BARE = details ex_qag ex_qags ex_triex ex_decuhr2d ex_cutet ex_decuhr3d simplexpapertest
DRIVERS_OBJ = $(DRIVERS_BARE:%=$(LOCAL_BUILD)/Drivers/%.o)
# details    : prints info on current distribution
# ex_qag     : 1-dimensional integration, tests of QAG from QUADPACK
# ex_qags    : 1-dimensional integration, tests of QAGS from QUADPACK
# ex_triex   : 2-dimensional integration (triangle), tests of TRIEX
# ex_decuhr2d: 2-dimensional integration (square), test of DECUHR
# ex_cutet   : 3-dimensional integration (tetrahedron), test of DCUTET
# ex_decuhr3d: 3-dimensional integration (cube), tests of DECUHR
# simplexpapertest : 5-dimension integration (simplex)


all: cubpack

cubpack: prebuild $(CORE_OBJ)
	$(AR) $(ARFLAGS) $(LIB)/libcubpack.a $(CORE_OBJ)
	$(RANLIB) $(LIB)/libcubpack.a
	
prebuild:
	cd $(LOCAL_BUILD); \
	$(MKDIR_P) Core Drivers

cubpack_drivers: cubpack $(DRIVERS_OBJ)
	$(foreach DRIVER,$(DRIVERS_BARE), $(FF) $(FCFLAGS) -o $(LOCAL_BUILD)/$(DRIVER) $(LOCAL_BUILD)/Drivers/$(DRIVER).o -L$(LIB) -lcubpack;)

$(LOCAL_BUILD)/Core/%.o: $(LOCAL_SRC)/Core/%.f90
	$(FC) -c $(FCFLAGS) -o $@ -module $(LOCAL_BUILD)/Core $<

$(LOCAL_BUILD)/Drivers/%.o: $(LOCAL_SRC)/Drivers/%.f90
	$(FC) -c $(FCFLAGS) -o $@ -module $(LOCAL_BUILD)/Drivers -I$(LOCAL_BUILD)/Core $<

clean:
	rm -f -r $(LOCAL_BUILD)/Core/*.o
	rm -f -r $(LOCAL_BUILD)/Core/*.mod
	rm -f -r $(LOCAL_BUILD)/Drivers/*.o
	rm -f -r $(LOCAL_BUILD)/Drivers/*.mod
	rm -f -r $(LOCAL_ROOT)/*.mod
	rm -f -r $(LOCAL_ROOT)/*_genmod*
	rm -f $(LIB)/libcubpack.a

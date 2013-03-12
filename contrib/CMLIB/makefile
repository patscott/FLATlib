# Makefile for CMLIB within flatlib
#
# Pat Scott Feb 2009
# pat@fysik.su.se

LOCAL_ROOT = $(PWD)
LOCAL_BUILD = ${LOCAL_ROOT}/build
LOCAL_SRC = ${LOCAL_ROOT}/src

BSPLINE_OBJ_BARE = dbintk.o dbnfac.o dbnslv.o dbspvn.o dbvalu.o dintrv.o
MACHCON_OBJ_BARE = i1mach.o
TENSBS_OBJ_BARE = db2ink.o db2val.o db3ink.o db3val.o dbknot.o dbtpcf.o
XERROR_OBJ_BARE = fdump.o j4save.o xerabt.o xerctl.o xerprt.o xerror.o xerrwv.o xersav.o xgetua.o

BSPLINE_OBJ = $(BSPLINE_OBJ_BARE:%.o=dbsplin/%.o)
MACHCON_OBJ = $(MACHCON_OBJ_BARE:%.o=machcon/%.o)
TENSBS_OBJ = $(TENSBS_OBJ_BARE:%.o=dtensbs/%.o)
XERROR_OBJ = $(XERROR_OBJ_BARE:%.o=xerror/%.o)

CMLIB_OBJ_BARE = $(BSPLINE_OBJ) $(MACHCON_OBJ) $(TENSBS_OBJ) $(XERROR_OBJ) 
CMLIB_OBJ = $(CMLIB_OBJ_BARE:%.o=$(LOCAL_BUILD)/%.o)


all: cmlib

cmlib:  prebuild $(CMLIB_OBJ)
	$(AR) $(ARFLAGS) $(LIB)/libcmlib.a $(CMLIB_OBJ)
	$(RANLIB) $(LIB)/libcmlib.a

prebuild:
	cd $(LOCAL_BUILD); \
	$(MKDIR_P) dbsplin machcon dtensbs xerror
	
$(LOCAL_BUILD)/%.o: $(LOCAL_SRC)/%.f
	$(FC) -c $(FCFLAGS) -o $@ $<

clean:
	rm -f -r $(LOCAL_BUILD)/bspline/*.o
	rm -f -r $(LOCAL_BUILD)/bspline/*.mod
	rm -f -r $(LOCAL_BUILD)/machcon/*.o
	rm -f -r $(LOCAL_BUILD)/machcon/*.mod
	rm -f -r $(LOCAL_BUILD)/tensbs/*.o
	rm -f -r $(LOCAL_BUILD)/tensbs/*.mod
	rm -f -r $(LOCAL_BUILD)/xerror/*.o
	rm -f -r $(LOCAL_BUILD)/xerror/*.mod
	rm -f -r $(LOCAL_ROOT)/*.mod
	rm -f -r $(LOCAL_ROOT)/*_genmod*
	rm -f $(LIB)/libcmlib.a

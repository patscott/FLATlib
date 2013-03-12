# Makefile for FITPACK within flatlib
#
# Pat Scott Mar 2009
# pat@fysik.su.se

LOCAL_ROOT = $(PWD)
LOCAL_BUILD = ${LOCAL_ROOT}/build
LOCAL_SRC = ${LOCAL_ROOT}/src

CORE_OBJ_BARE = fitpack.o
CORE_OBJ = $(CORE_OBJ_BARE:%.o=$(LOCAL_BUILD)/%.o)

all: fitpack

fitpack: $(CORE_OBJ)
	$(AR) $(ARFLAGS) $(LIB)/libfitpack.a $(CORE_OBJ)
	$(RANLIB) $(LIB)/libfitpack.a
	
$(LOCAL_BUILD)/%.o: $(LOCAL_SRC)/%.f
	$(FC) -c $(FCFLAGS) -o $@ $<

clean:
	rm -f -r $(LOCAL_BUILD)/*.o
	rm -f -r $(LOCAL_BUILD)/*.mod
	rm -f -r $(LOCAL_ROOT)/*.mod
	rm -f -r $(LOCAL_ROOT)/*_genmod*
	rm -f $(LIB)/libfitpack.a

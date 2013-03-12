# Makefile for NSWC library within flatlib
#
# Pat Scott Mar 2009
# pat@fysik.su.se

LOCAL_ROOT = $(PWD)
LOCAL_BUILD = ${LOCAL_ROOT}/build
LOCAL_SRC = ${LOCAL_ROOT}/src

CORE_OBJ_BARE = erfinv.o
CORE_OBJ = $(CORE_OBJ_BARE:%.o=$(LOCAL_BUILD)/%.o)

all: nswc

nswc: $(CORE_OBJ)
	$(AR) $(ARFLAGS) $(LIB)/libnswc.a $(CORE_OBJ)
	$(RANLIB) $(LIB)/libnswc.a
	
$(LOCAL_BUILD)/%.o: $(LOCAL_SRC)/%.f90
	$(FC) -c $(FCFLAGS) -o $@ $<

clean:
	rm -f -r $(LOCAL_BUILD)/*.o
	rm -f -r $(LOCAL_BUILD)/*.mod
	rm -f -r $(LOCAL_ROOT)/*.mod
	rm -f -r $(LOCAL_ROOT)/*_genmod*
	rm -f $(LIB)/libnswc.a

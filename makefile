# Makefile for flatlib v0.6
#
# Note that flatlib requires
# the cfitsio, fftw3 and m
# libraries to be installed,
# and the header fftw3.f must
# be in an included path.
#
#
# Pat Scott Oct 2012
# pat@fysik.su.se

export

FOPT=-O3 -extend_source -parallel -openmp -lpthread
# non-MPI mode
FC=ifort
FCFLAGS=$(FOPT) -warn #-check all
# MPI mode - not implemented yet in the code
#FF = mpif90
#FFLAGS=$(FOPT) -warn -DMPI #-check all
RANLIB = ranlib
MKDIR_P = mkdir -p
AR = ar
ARFLAGS = rv

FLAT_ROOT = $(PWD)
LIB = ${FLAT_ROOT}/lib
SRC = ${FLAT_ROOT}/src
BUILD = ${FLAT_ROOT}/build
DATA = ${FLAT_ROOT}/data
CONTRIB = ${FLAT_ROOT}/contrib

CUBPACKDIR = $(CONTRIB)/CUBPACK
FITPACKDIR = $(CONTRIB)/FITPACK
CMLIBDIR = $(CONTRIB)/CMLIB
NSWCDIR = $(CONTRIB)/NSWC

COMMONOBJ_BARE = flatCommon.o

OBJ_BARE = flatUtils.o flatIRFparams.o flatIRFs.o\
           flatFFTW.o flatPrecompute.o flatIRFini.o\
           flatConvolve_fast.o

TESTOBJ_BARE = flatTest.o
AVEROBJ_BARE = flatAverage.o

COMMONOBJ = $(COMMONOBJ_BARE:%.o=$(BUILD)/%.o)
OBJ = $(OBJ_BARE:%.o=$(BUILD)/%.o)
TESTOBJ = $(TESTOBJ_BARE:%.o=$(BUILD)/%.o)
AVEROBJ = $(AVEROBJ_BARE:%.o=$(BUILD)/%.o)


all: flattest flataverage

flatlib: cmlib cubpack fitpack nswc common $(OBJ) 
	$(AR) $(ARFLAGS) $(LIB)/libflatlib.a $(COMMONOBJ) $(OBJ)
	cd $(LIB); $(AR) x libcmlib.a; $(AR) x libcubpack.a; $(AR) x libfitpack.a; $(AR) x libnswc.a
	$(AR) $(ARFLAGS) $(LIB)/libflatlib.a $(LIB)/*.o
	rm -f $(LIB)/*.o
	$(RANLIB) $(LIB)/libflatlib.a

cmlib:
	cd $(CMLIBDIR); $(MAKE)
	
cubpack:
	cd $(CUBPACKDIR); $(MAKE)

cubpack_drivers:
	cd $(CUBPACKDIR); $(MAKE) cubpack_drivers

fitpack:
	cd $(FITPACKDIR); $(MAKE) fitpack

nswc:
	cd $(NSWCDIR); $(MAKE) nswc

flattest: flatlib $(TESTOBJ)
	$(FC) $(FCFLAGS) -L$(LIB) -I$(BUILD) -o $(FLAT_ROOT)/flattest \
	$(TESTOBJ) -lflatlib -lcfitsio -lfftw3 -lm

flataverage: flatlib $(AVEROBJ)
	$(FC) $(FCFLAGS) -L$(LIB) -I$(BUILD) -o $(DATA)/flataverage \
	$(AVEROBJ) -lflatlib -lcfitsio -lfftw3 -lm

common: $(COMMONOBJ_BARE:%.o=$(SRC)/%.f90)
	cd $(BUILD); \
	$(FC) -c $(FCFLAGS) -I$(CONTRIB)/CUBPACK/build/Core $?

$(SRC)/flatCommon.f90: $(SRC)/flatCommon.f90.template
	perl $(FLAT_ROOT)/scr/config.pl $(FLAT_ROOT)

$(BUILD)/%.o: $(SRC)/%.f90 $(COMMONOBJ_BARE:%.o=$(SRC)/%.f90)
	cd $(BUILD); \
	$(FC) -c -o $@ $(FCFLAGS) -I$(CONTRIB)/CUBPACK/build/Core $<

$(BUILD)/%.o: $(SRC)/%.f $(COMMONOBJ_BARE:%.o=$(SRC)/%.f90)
	cd $(BUILD); \
	$(FC) -c -o $@ $(FCFLAGS) -I$(CONTRIB)/CUBPACK/build/Core $<

clean:
	rm -f $(BUILD)/*
	rm -f $(LIB)/*
	rm -f -r $(CONTRIB)/*/build/* 
	rm -f $(DATA)/flataverage
	rm -f flattest
	rm -f $(SRC)/flatCommon.f90
	cd $(CUBPACKDIR); $(MAKE) clean
	cd $(CMLIBDIR); $(MAKE) clean
	cd $(FITPACKDIR); $(MAKE) clean
	cd $(NSWCDIR); $(MAKE) clean

cleanall: clean
	rm -f $(DATA)/*/*_mean.dat


FORTRAN=gnu
DEBUG=false
## intel 
#INCLUDE=-I/opt/intel_compiled/include
#LIB_INC=-L/opt/intel_compiled/lib
CPU=native
INCLUDE=-I/usr/local/include
LIB_INC=-L/usr/local/lib
F90SAC = f90sac.o f90sac_csubs.o

# --------------------- don't touch this -------------------------------#

FFTLIB=-lfftw3 -lm
NCLIB=-lnetcdff -lnetcdf

ifeq ($(FORTRAN),gnu)
        FC=gfortran
        CC=gcc-11
        ifeq ($(DEBUG),true)
        	FFLAGS =-g -fbounds-check -Wall -fbacktrace -finit-real=nan $(INCLUDE)
        else
        	FFLAGS=-m64 -O3 -march=$(CPU) $(INCLUDE) -fopenmp
        endif
endif

ifeq ($(FORTRAN),intel)
        FC=ifort
        CC=icc
        ifeq ($(DEBUG),true)
		FFLAGS=-g -traceback -O0 -check bounds -warn all -debug extended -save-temps $(INCLUDE)
        else
        	FFLAGS=-m64 -march=$(CPU) $(INCLUDE) -qopenmp
        	CFLAGS=-m64 -fast -march=$(CPU)
        endif
endif

ifeq ($(FORTRAN),portland)
        FC=pgf95
        CC=pgcc
        ifeq ($(DEBUG),true)
        	FFLAGS=-g -m64 -traceback -Mbounds $(INCLUDE)
        else
        	FFLAGS=-m64 -O3 -fastsse $(INCLUDE) -mp
        endif
endif


%.o : %.c
	$(CC) -c $(CFLAGS) $<

%.o : %.f90
	$(FC) -c $(FFLAGS) $<

%.o : %.F90
	$(FC) -c $(FFLAGS) $<

OBJECTS = string_utility.o sa_io.o sa_math.o sa_fourier.o sa_array.o

all: tfreqfisher freqfisher timefisher timefisherCDF tdoa tdoa_test ccts

ccts: ${F90SAC} ${OBJECTS} ccts.o
		$(FC) ${F90SAC} ${OBJECTS} ccts.o \
		$(FFLAGS) -o ccts $(LIB_INC) $(FFTLIB)
#		install_name_tool -change @rpath/libiomp5.dylib $(OMPPATH)/libiomp5.dylib ccts

tdoa_test: ${OBJECTS} tdoa_test.o
		$(FC) sa_math.o tdoa_test.o \
		$(FFLAGS) -o tdoa_test

tdoa: ${F90SAC} ${OBJECTS} tdoa.o
		$(FC) ${F90SAC} ${OBJECTS} tdoa.o \
		$(FFLAGS) -o tdoa $(LIB_INC) $(FFTLIB)
#		install_name_tool -change @rpath/libiomp5.dylib $(OMPPATH)/libiomp5.dylib tdoa

timefisher: ${F90SAC} ${OBJECTS} timefisher.o
		$(FC) ${F90SAC} ${OBJECTS} timefisher.o \
		$(FFLAGS) -o timefisher $(LIB_INC) $(FFTLIB)
#		install_name_tool -change @rpath/libiomp5.dylib $(OMPPATH)/libiomp5.dylib timefisher

timefisherCDF: ${F90SAC} ${OBJECTS} timefisherCDF.o
		$(FC) ${F90SAC} ${OBJECTS} timefisherCDF.o \
		$(FFLAGS) -o timefisherCDF $(LIB_INC) $(FFTLIB) $(NCLIB)
		#install_name_tool -change @rpath/libiomp5.dylib $(OMPPATH)/libiomp5.dylib timefisherCDF

freqfisher: ${F90SAC} ${OBJECTS} freqfisher.o
		$(FC) ${F90SAC} ${OBJECTS} freqfisher.o \
		$(FFLAGS) -o freqfisher $(LIB_INC) $(FFTLIB)
#		install_name_tool -change @rpath/libiomp5.dylib $(OMPPATH)/libiomp5.dylib freqfisher

tfreqfisher: ${F90SAC} ${OBJECTS} tfreqfisher.o
		$(FC) ${F90SAC} ${OBJECTS} tfreqfisher.o \
		$(FFLAGS) -o tfreqfisher $(LIB_INC) $(FFTLIB)
#		install_name_tool -change @rpath/libiomp5.dylib $(OMPPATH)/libiomp5.dylib tfreqfisher

clean:
	rm -rf *.o *.mod *.i90 tfreqfisher freqfisher timefisher timefisherCDF tdoa tdoa_test ccts

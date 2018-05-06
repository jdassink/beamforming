FORTRAN=gnu
DEBUG=false
INCLUDE=-I/opt/local/include
LIB_INC=-L/opt/local/lib
CPU=corei7
#CPU=native
#INCLUDE=-I/usr/include
#LIB_INC=-L/usr/lib64
F90SAC = f90sac.o f90sac_csubs.o

# --------------------- don't touch this -------------------------------#

FFTLIB=-lfftw3 -lm

ifeq ($(FORTRAN),gnu)
        FC=gfortran
        CC=gcc
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
        	FFLAGS=-m64 -fast -march=$(CPU) $(INCLUDE) -qopenmp
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

all: freqfisher timefisher tdoa ccts

ccts: ${F90SAC} ${OBJECTS} ccts.o
		$(FC) ${F90SAC} ${OBJECTS} ccts.o \
		$(FFLAGS) -o ccts $(LIB_INC) $(FFTLIB)

tdoa: ${F90SAC} ${OBJECTS} tdoa.o
		$(FC) ${F90SAC} ${OBJECTS} tdoa.o \
		$(FFLAGS) -o tdoa $(LIB_INC) $(FFTLIB)

timefisher: ${F90SAC} ${OBJECTS} timefisher.o
		$(FC) ${F90SAC} ${OBJECTS} timefisher.o \
		$(FFLAGS) -o timefisher $(LIB_INC) $(FFTLIB)

freqfisher: ${F90SAC} ${OBJECTS} freqfisher.o
		$(FC) ${F90SAC} ${OBJECTS} freqfisher.o \
		$(FFLAGS) -o freqfisher $(LIB_INC) $(FFTLIB)

clean:
	rm -rf *.o *.mod *.i90 freqfisher timefisher tdoa ccts

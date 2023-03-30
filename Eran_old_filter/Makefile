SHELL = /bin/sh
OS=Linux

Linux_LIB = -shared-intel -Wl,--start-group /opt/intel/oneapi/mkl/latest/lib/intel64/libmkl_intel_lp64.a /opt/intel/oneapi/mkl/latest/lib/intel64/libmkl_sequential.a /opt/intel/oneapi/mkl/latest/lib/intel64/libmkl_core.a -lfftw3_threads -lfftw3 -lpthread
#Linux_LIB = -llapack -lblas -lfftw3 -lm
Linux_MYLIB = 
Linux_CLIB = -lm
Linux_FCLIB = 	

Linux_DYLIB = -shared-intel -Wl,--start-group /opt/intel/oneapi/mkl/latest/lib/intel64/libmkl_intel_lp64.a /opt/intel/oneapi/mkl/latest/lib/intel64/libmkl_sequential.a /opt/intel/oneapi/mkl/latest/lib/intel64/libmkl_core.a -lfftw3_threads -lfftw3 -lpthread
CYGWIN_NT-6.3_LIB = -llapack -lblas -lfftw3_threads -lfftw3 -lpthread
CYGWIN_NT-6.3_MYLIB = 
CYGWIN_NT-6.3_CLIB =  -lm
CYGWIN_NT-6.3_FCLIB = 	

LIB = $(${OS}_DYLIB) $(${OS}_CLIB) $(${OS}_FCLIB) $(${OS}_MYLIB)

# flags ...
#Linux_OPTFLAGS = -DMKL_ILP64 -fast -qopenmp -static
#Linux_OPTFLAGS = -DMKL_ILP64 -fast -qopenmp
Linux_OPTFLAGS = -fast -qopenmp -w

#Linux_OPTFLAGS = -DMKL_ILP64 -O2 -xSSE4.1 -qopenmp

Linux_CFLAGS = $(${OS}_OPTFLAGS) 
Linux_FFLAGS = $(${OS}_OPTFLAGS)

CYGWIN_NT-6.3_OPTFLAGS = -O3 -ffast-math 
CYGWIN_NT-6.3_CFLAGS = $(${OS}_OPTFLAGS) -DFFT_FFTW
CYGWIN_NT-6.3_FFLAGS = $(${OS}_OPTFLAGS) 

MAINNAM = filter

#compiler ...
Linux_CC = icx
Linux_FF = ifx
Linux_LD = icx
CYGWIN_NT-6.3_CC = gcc
CYGWIN_NT-6.3_FF = gfortran
CYGWIN_NT-6.3_LD = gcc

OBJECTS = \
	main.o init.o size.o norm.o nerror.o coeff.o read.o interpolate.o states.o \
	ortho.o energy.o hamiltonian.o hnorm.o filter.o Hmat.o absorption.o rand.o


# compilation ...

.f.o:
	$(${OS}_FF) $(${OS}_FFLAGS) -c  $*.f
.c.o:
	$(${OS}_CC) -DOS_$(OS) $(${OS}_CFLAGS) -c  $*.c  

$(MAINNAM): $(OBJECTS) 
	$(${OS}_LD) -o $(MAINNAM).x $(${OS}_CFLAGS) $(OBJECTS) $(LIB)

clean:
	/bin/rm *.o *.x

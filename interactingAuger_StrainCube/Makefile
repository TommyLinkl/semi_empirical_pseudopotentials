SHELL = /bin/sh

Linux_LIB = -shared-intel -Wl,--start-group /opt/intel/compilers_and_libraries_2018/linux/mkl/lib/intel64/libmkl_intel_lp64.a /opt/intel/compilers_and_libraries_2018/linux/mkl/lib/intel64/libmkl_sequential.a /opt/intel/compilers_and_libraries_2018/linux/mkl/lib/intel64/libmkl_core.a -lfftw3 -fopenmp
Linux_LIB = -shared-intel -Wl,--start-group /opt/intel/compilers_and_libraries_2018/linux/mkl/lib/intel64/libmkl_intel_ilp64.a /opt/intel/compilers_and_libraries_2018/linux/mkl/lib/intel64/libmkl_lapack95_ilp64.a /opt/intel/compilers_and_libraries_2018/linux/mkl/lib/intel64/libmkl_blas95_ilp64.a /opt/intel/compilers_and_libraries_2018/linux/mkl/lib/intel64/libmkl_intel_thread.a /opt/intel/compilers_and_libraries_2018/linux/mkl/lib/intel64/libmkl_core.a -Wl,--end-group -qopenmp -lpthread -lfftw3 
Linux_LIB = -shared-intel -Wl,--start-group /opt/intel/mkl/lib/intel64/libmkl_intel_ilp64.a /opt/intel/mkl/lib/intel64/libmkl_lapack95_ilp64.a /opt/intel/mkl/lib/intel64/libmkl_blas95_ilp64.a /opt/intel/mkl/lib/intel64/libmkl_intel_thread.a /opt/intel/mkl/lib/intel64/libmkl_core.a -Wl,--end-group -qopenmp -lpthread -lfftw3 
Linux_MYLIB = 
Linux_CLIB = -lm
Linux_FCLIB = 	

CYGWIN_NT-10.0_LIB = -lblas -llapack -lfftw3_threads -lfftw3 -lpthread
CYGWIN_NT-10.0_MYLIB = 
CYGWIN_NT-10.0_CLIB =  -lm
CYGWIN_NT-10.0_FCLIB = 	

LIB = $(${OS}_LIB) $(${OS}_CLIB) $(${OS}_FCLIB) $(${OS}_MYLIB)

# flags ...
# Linux_OPTFLAGS = -fast -qopenmp
# Linux_OPTFLAGS = -DMKL_ILP64 -O3 -xSSE3 -no-prec-div -fp-model fast=2 -qopenmp -static -g
Linux_OPTFLAGS = -DMKL_ILP64 -O2 -xSSE4.1 -qopenmp -g
Linux_CFLAGS = $(${OS}_OPTFLAGS) -DFFT_FFTW
Linux_FFLAGS = $(${OS}_OPTFLAGS)

CYGWIN_NT-10.0_OPTFLAGS = -O3 -ffast-math
CYGWIN_NT-10.0_CFLAGS = $(${OS}_OPTFLAGS) -DFFT_FFTW
CYGWIN_NT-10.0_FFLAGS = $(${OS}_OPTFLAGS) 

MAINNAM = auger-bse-boltz

#compiler ...
Linux_CC = icx
Linux_FF = ifx
Linux_LD = icx

CYGWIN_NT-10.0_CC = gcc
CYGWIN_NT-10.0_FF = gfortran
CYGWIN_NT-10.0_LD = gcc

OBJECTS = \
	main.o vector.o init.o write.o size.o norm.o nerror.o read.o gauss.o hartree.o auger.o energy.o interpolate.o hamiltonian.o hnorm.o filter.o coeff.o rand.o ortho.o Hmat.o exciton.o generate-filter-states.o coulomb.o states.o


# compilation ...

.f.o:
	$(${OS}_FF) $(${OS}_FFLAGS) -c  $*.f
.c.o:
	$(${OS}_CC) -DOS_$(OS) $(${OS}_CFLAGS) -c  $*.c  

$(MAINNAM): $(OBJECTS) 
	$(${OS}_LD) -o $(MAINNAM).x $(${OS}_CFLAGS) $(OBJECTS) $(LIB)

clean:
	/bin/rm *.o *.x

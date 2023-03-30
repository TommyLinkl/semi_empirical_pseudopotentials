SHELL = /bin/sh

Linux_LIB = -shared-intel -Wl,--start-group /opt/intel/compilers_and_libraries_2018/linux/mkl/lib/intel64/libmkl_intel_lp64.a /opt/intel/compilers_and_libraries_2018/linux/mkl/lib/intel64/libmkl_sequential.a /opt/intel/compilers_and_libraries_2018/linux/mkl/lib/intel64/libmkl_core.a -lfftw3 -fopenmp
Linux_LIB = -shared-intel -Wl,--start-group /opt/intel/compilers_and_libraries_2018/linux/mkl/lib/intel64/libmkl_intel_ilp64.a /opt/intel/compilers_and_libraries_2018/linux/mkl/lib/intel64/libmkl_lapack95_ilp64.a /opt/intel/compilers_and_libraries_2018/linux/mkl/lib/intel64/libmkl_blas95_ilp64.a /opt/intel/compilers_and_libraries_2018/linux/mkl/lib/intel64/libmkl_intel_thread.a /opt/intel/compilers_and_libraries_2018/linux/mkl/lib/intel64/libmkl_core.a -Wl,--end-group -qopenmp -lpthread -lfftw3
Linux_MYLIB = 
Linux_CLIB = -lm
Linux_FCLIB = 	

CYGWIN_NT-6.3_LIB = -llapack -lblas -lfftw3_threads -lfftw3 -lpthread
CYGWIN_NT-6.3_MYLIB = 
CYGWIN_NT-6.3_CLIB =  -lm
CYGWIN_NT-6.3_FCLIB = 	

LIB = $(${OS}_LIB) $(${OS}_CLIB) $(${OS}_FCLIB) $(${OS}_MYLIB)

# flags ...
#Linux_OPTFLAGS = -fast -qopenmp -static
Linux_OPTFLAGS = -DMKL_ILP64 -O3 -xSSE3 -no-prec-div -fp-model fast=2 -qopenmp -static
#Linux_OPTFLAGS = -O3 -static
Linux_CFLAGS = $(${OS}_OPTFLAGS) -DFFT_FFTW -I /opt/rocks/include
Linux_FFLAGS = $(${OS}_OPTFLAGS)

CYGWIN_NT-6.3_OPTFLAGS = -O3 -ffast-math 
CYGWIN_NT-6.3_CFLAGS = $(${OS}_OPTFLAGS) -DFFT_FFTW
CYGWIN_NT-6.3_FFLAGS = $(${OS}_OPTFLAGS) 

MAINNAM = spinOrbitFilter

#compiler ...
Linux_CC = icc
Linux_FF = ifort
Linux_LD = icc
CYGWIN_NT-6.3_CC = gcc
CYGWIN_NT-6.3_FF = gfortran
CYGWIN_NT-6.3_LD = gcc

OBJECTS = \
	main.o init.o size.o norm.o nerror.o coeff.o read.o interpolate.o \
	ortho.o energy.o hamiltonian.o hnorm.o filter.o Hmat.o dipole.o rand.o \
	num_int.o 


# compilation ...

.f.o:
	$(${OS}_FF) $(${OS}_FFLAGS) -c  $*.f
.c.o:
	$(${OS}_CC) -DOS_$(OS) $(${OS}_CFLAGS) -c  $*.c  

$(MAINNAM): $(OBJECTS) 
	$(${OS}_LD) -o $(MAINNAM).x $(${OS}_CFLAGS) $(OBJECTS) $(LIB)

clean:
	/bin/rm *.o *.x

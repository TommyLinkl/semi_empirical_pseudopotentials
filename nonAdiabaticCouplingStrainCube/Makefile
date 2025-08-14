SHELL = /bin/sh

# libraries
Linux_LIB = -shared-intel -Wl,--start-group /opt/intel/compilers_and_libraries_2018/linux/mkl/lib/intel64/libmkl_intel_lp64.a /opt/intel/compilers_and_libraries_2018/linux/mkl/lib/intel64/libmkl_sequential.a /opt/intel/compilers_and_libraries_2018/linux/mkl/lib/intel64/libmkl_core.a -lfftw3 -fopenmp
Linux_LIB = -shared-intel -Wl,--start-group /opt/intel/compilers_and_libraries_2018/linux/mkl/lib/intel64/libmkl_intel_ilp64.a /opt/intel/compilers_and_libraries_2018/linux/mkl/lib/intel64/libmkl_lapack95_ilp64.a /opt/intel/compilers_and_libraries_2018/linux/mkl/lib/intel64/libmkl_blas95_ilp64.a /opt/intel/compilers_and_libraries_2018/linux/mkl/lib/intel64/libmkl_intel_thread.a /opt/intel/compilers_and_libraries_2018/linux/mkl/lib/intel64/libmkl_core.a -Wl,--end-group -qopenmp -lpthread -lfftw3 
Linux_MYLIB = 
Linux_CLIB = -lm
Linux_FCLIB = 	

LIB = $(${OS}_LIB) $(${OS}_CLIB) $(${OS}_FCLIB) $(${OS}_MYLIB)

# flags 
Linux_OPTFLAGS = -fast -qopenmp
Linux_OPTFLAGS = -DMKL_ILP64 -O3 -xSSE3 -no-prec-div -fp-model fast=2 -qopenmp # -static
Linux_CFLAGS = $(${OS}_OPTFLAGS) -DFFT_FFTW
Linux_FFLAGS = $(${OS}_OPTFLAGS)

MAINNAM = nonadiabaticStrainCube

# compiler 
Linux_CC = icc
Linux_FF = ifort
Linux_LD = icc

OBJECTS = \
	main.o read.o write.o grid.o size.o vector.o errorHandling.o pseudopotential.o nonadiabatic.o 

# compilation ...
.f.o:
	$(${OS}_FF) $(${OS}_FFLAGS) -c  $*.f
.c.o:
	$(${OS}_CC) -DOS_$(OS) $(${OS}_CFLAGS) -c  $*.c  

$(MAINNAM): $(OBJECTS) 
	$(${OS}_LD) -o $(MAINNAM).x $(${OS}_CFLAGS) $(OBJECTS) $(LIB)

clean:
	/bin/rm *.o *.x

SHELL = /bin/sh

# libraries ...
Linux_LIB = -shared-intel -Wl,--start-group /opt/intel/compilers_and_libraries_2018/linux/mkl/lib/intel64/libmkl_intel_lp64.a /opt/intel/compilers_and_libraries_2018/linux/mkl/lib/intel64/libmkl_sequential.a /opt/intel/compilers_and_libraries_2018/linux/mkl/lib/intel64/libmkl_core.a -lfftw3 -fopenmp
Linux_LIB = -shared-intel -Wl,--start-group /opt/intel/compilers_and_libraries_2018/linux/mkl/lib/intel64/libmkl_intel_ilp64.a /opt/intel/compilers_and_libraries_2018/linux/mkl/lib/intel64/libmkl_lapack95_ilp64.a /opt/intel/compilers_and_libraries_2018/linux/mkl/lib/intel64/libmkl_blas95_ilp64.a /opt/intel/compilers_and_libraries_2018/linux/mkl/lib/intel64/libmkl_intel_thread.a /opt/intel/compilers_and_libraries_2018/linux/mkl/lib/intel64/libmkl_core.a -Wl,--end-group -qopenmp -lpthread -lfftw3
Linux_CLIB = -lm
Linux_DYLIB = -lmkl_intel_ilp64 -lmkl_lapack95_ilp64 -lmkl_blas95_ilp64 -lmkl_intel_thread -lmkl_core -lm

LIB = $(${OS}_DYLIB) $(${OS}_LIB) $(${OS}_CLIB) $(${OS}_FCLIB) $(${OS}_MYLIB)

# flags ...
Linux_OPTFLAGS = -DMKL_ILP64 -O2 -xSSE4.1 
Linux_CFLAGS = $(${OS}_OPTFLAGS) -qopenmp
Linux_FFLAGS = $(${OS}_OPTFLAGS)

MAINNAM = filter

# compiler ...
Linux_CC = icx
Linux_FF = ifort
Linux_LD = icx

OBJECTS = \
	main.o write.o init.o size.o norm.o nerror.o coeff.o read.o interpolate.o \
	ortho.o energy.o hamiltonian.o hnorm.o filter.o Hmat.o dipole.o rand.o

# compilation ...
.f.o:
	$(${OS}_FF) $(${OS}_FFLAGS) -c  $*.f
.c.o:
	$(${OS}_CC) -DOS_$(OS) $(${OS}_CFLAGS) -c  $*.c  

$(MAINNAM): $(OBJECTS) 
	$(${OS}_LD) -o $(MAINNAM).x $(${OS}_CFLAGS) $(OBJECTS) $(LIB)

clean:
	/bin/rm *.o *.x

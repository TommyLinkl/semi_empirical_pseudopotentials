SHELL = /bin/sh
OS = Linux

# libraries ...
Linux_LIB = -shared-intel -Wl,--start-group /opt/intel/oneapi/mkl/latest/lib/intel64/libmkl_intel_lp64.a /opt/intel/oneapi/mkl/latest/lib/intel64/libmkl_sequential.a /opt/intel/oneapi/mkl/latest/lib/intel64/libmkl_core.a -lfftw3_threads -lfftw3 -lpthread
Linux_CLIB = -lm
Linux_DYLIB = -shared-intel -Wl,--start-group /opt/intel/oneapi/mkl/latest/lib/intel64/libmkl_intel_lp64.a /opt/intel/oneapi/mkl/latest/lib/intel64/libmkl_sequential.a /opt/intel/oneapi/mkl/latest/lib/intel64/libmkl_core.a -lfftw3_threads -lfftw3 -lpthread

LIB = $(${OS}_DYLIB) $(${OS}_LIB) $(${OS}_CLIB) $(${OS}_FCLIB) $(${OS}_MYLIB)

# flags 
Linux_OPTFLAGS = -fast
Linux_CFLAGS = $(${OS}_OPTFLAGS) -qopenmp
Linux_FFLAGS = $(${OS}_OPTFLAGS)

MAINNAM = filter_debug-strain

# compiler ...
Linux_CC = icx
Linux_FF = ifx
Linux_LD = icx

OBJECTS = \
	main.o vector.o write.o init.o size.o norm.o nerror.o coeff.o read.o interpolate.o \
	ortho.o energy.o hamiltonian.o hnorm.o filter.o Hmat.o dipole.o rand.o slice_psi.o

# compilation ...
.f.o:
	$(${OS}_FF) $(${OS}_FFLAGS) -c  $*.f
.c.o:
	$(${OS}_CC) -DOS_$(OS) $(${OS}_CFLAGS) -c  $*.c  

$(MAINNAM): $(OBJECTS) 
	$(${OS}_LD) -o $(MAINNAM).x $(${OS}_CFLAGS) $(OBJECTS) $(LIB)

clean:
	/bin/rm *.o *.x

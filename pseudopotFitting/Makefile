# Source code 
OBJS = main.o init.o pot.o read.o states.o hamiltonian.o diag.o write.o vector.o fit.o rand.o

# Tells the program to use the intel c compiler - faster than the gnu compiler 
CC = icc
FF = ifort
LD = icc

# Gives the path to the intel compilers needed in this program
LIB = -lm 
LIB = -shared-intel -Wl,--start-group /opt/intel/compilers_and_libraries_2018/linux/mkl/lib/intel64/libmkl_intel_lp64.a /opt/intel/compilers_and_libraries_2018/linux/mkl/lib/intel64/libmkl_sequential.a /opt/intel/compilers_and_libraries_2018/linux/mkl/lib/intel64/libmkl_core.a -Wl,--end-group -qopenmp -lpthread -lfftw3
LIB = -shared-intel -Wl,--start-group /opt/intel/compilers_and_libraries_2018/linux/mkl/lib/intel64/libmkl_intel_ilp64.a /opt/intel/compilers_and_libraries_2018/linux/mkl/lib/intel64/libmkl_lapack95_ilp64.a /opt/intel/compilers_and_libraries_2018/linux/mkl/lib/intel64/libmkl_blas95_ilp64.a /opt/intel/compilers_and_libraries_2018/linux/mkl/lib/intel64/libmkl_intel_thread.a /opt/intel/compilers_and_libraries_2018/linux/mkl/lib/intel64/libmkl_core.a -Wl,--end-group -qopenmp -lpthread -lfftw3 -mkl 
OPTFLAGS = -DMKL_ILP64 -O2 -xSSE3 -no-prec-div -fp-model fast=2 -qopenmp -static 
CFLAGS = $(OPTFLAGS) -DFFT_FFTW

fit: $(OBJS)
	$(CC) $(CFLAGS) -o fit.x $(OBJS) $(LIB)
clean:
	/bin/rm *.o *.x


OBJS = main.o init.o read.o pot.o states.o hamiltonian.o diag.o write.o vector.o fit.o rand.o

# Tells the program to use the intel c compiler - faster than the gnu compiler 
CC = icc

# Gives the path to the intel compilers needed in this program
LIB = -shared-intel -Wl,--start-group /opt/intel/compilers_and_libraries_2018/linux/mkl/lib/intel64/libmkl_intel_lp64.a /opt/intel/compilers_and_libraries_2018/linux/mkl/lib/intel64/libmkl_sequential.a /opt/intel/compilers_and_libraries_2018/linux/mkl/lib/intel64/libmkl_core.a -Wl,--end-group -qopenmp -lpthread -lfftw3 
CFLAGS = -fast -qopenmp  

fit: $(OBJS)
	$(CC) $(CFLAGS) -o soFit.x $(OBJS) $(LIB)
clean:
	/bin/rm *.o *.x


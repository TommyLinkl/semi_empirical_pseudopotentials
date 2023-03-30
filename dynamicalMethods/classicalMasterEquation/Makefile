OBJS = main.o traj.o read.o rand.o write.o init.c

# Tells the program to use the intel c compiler - faster than the gnu compiler 
CC = icc

# Gives the path to the intel compilers needed in this program
#LIB = -shared-intel -Wl,--start-group /opt/intel/compilers_and_libraries_2018/linux/mkl/lib/intel64/libmkl_intel_lp64.a /opt/        intel/compilers_and_libraries_2018/linux/mkl/lib/intel64/libmkl_sequential.a /opt/intel/compilers_and_libraries_2018/linux/mkl/lib/        intel64/libmkl_core.a -lfftw3 -fopenmp
#LIB = -shared-intel -Wl,--start-group /opt/intel/compilers_and_libraries_2018/linux/mkl/lib/intel64/libmkl_intel_ilp64.a /opt/intel/compilers_and_libraries_2018/linux/mkl/lib/intel64/libmkl_lapack95_ilp64.a /opt/intel/compilers_and_libraries_2018/linux/mkl/lib/intel64/libmkl_blas95_ilp64.a /opt/intel/compilers_and_libraries_2018/linux/mkl/lib/intel64/libmkl_intel_thread.a /opt/intel/compilers_and_libraries_2018/linux/mkl/lib/intel64/libmkl_core.a -Wl,--end-group -qopenmp -lpthread -lfftw3
CFLAGS = -O3 -fast #-qopenmp 

fit: $(OBJS)
	$(CC) $(CFLAGS) -o kinMCvDt.x $(OBJS) $(LIB)
clean:
	/bin/rm *.o *.x



CPP = g++
F90 = gfortran
F2P = python3 -m numpy.f2py

c_sources = source/lapjv.cpp
f_sources = source/lap.f90 source/utils.f90 source/tiling.f90 source/potential.f90 source/transform.f90 source/lapjv.f90
functions = munkres free_trans rot_mat center eye norm split det sort sphere circle distance derivative closest fixed_tmat fixed_tmat_int intoptimization fastoptimization lapjv optimize_vec
flags = --f90exec=gfortran --f90flags="-g -Og -fbacktrace -fopenmp" -lgomp 

#-I/usr/lib64/openmpi/lib/ -L/usr/lib64/openmpi/lib/ -lmpi

all: libs fmodules

fmodules: $(files)
	$(F2P) -c utils.o lapjv.o $(flags) -m p2ptrans.$@ $(files) only: $(functions)

libs: $(files) source/utils.f90 source/lapjv.cpp
	$(F90) -c -fPIC -fopenmp -Wall -g -Og source/utils.f90
	$(CPP) -fPIC -Wall -g -Og -c source/lapjv.cpp
clean:
	rm -f *.mod *.so *.a *.o

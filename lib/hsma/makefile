# compiler, and linking from C, fortran
CC=gcc # or icc if Intel package is used
CXX=g++
FC=gfortran

# set compiler flags for c and fortran
FFLAGS= -fPIC -O3 -march=native -funroll-loops -std=legacy -fopenmp 
CFLAGS= -fPIC -O3 -march=native -funroll-loops -std=c99 
CXXFLAGS= -std=c++11 -DSCTL_PROFILE=-1 -fPIC -O3 -march=native -funroll-loops

# set linking libraries
CLIBS = -lgfortran -lm -ldl 
LIBS = -lm -lstdc++ -fopenmp


# extra flags for multithreaded: C/Fortran, MATLAB
OMPFLAGS =-fopenmp
OMPLIBS =-lgomp 

# Test objects
TOBJS = ./hkrand.o ./dlaran.o
COMOBJS = ./besseljs3d.o ./cdjseval3d.o ./dfft.o \
	./fmmcommon.o ./legeexps.o ./prini.o \
	./rotgen.o ./rotproj.o ./rotviarecur.o \
	./tree_routs3d.o ./pts_tree3d.o ./yrecursion.o \
	./cumsum.o

LOBJS = ./lwtsexp_sep1.o ./l3dterms.o ./l3dtrans.o \
	./laprouts3d.o ./lfmm3d.o ./lfmm3dwrap.o \
	./lfmm3dwrap_legacy.o ./lfmm3dwrap_vec.o ./lwtsexp_sep2.o \
	./lpwrouts.o ./lapkernels.o ./lndiv.o ./cprini.o ./utils.o

objects := $(patsubst %.f,%.o,$(wildcard *.f))

objectsc := $(patsubst %.c,%.o,$(wildcard *.c))

# implicit rules for objects (note -o ensures writes to correct dir)
%.o: %.cpp %.h
	$(CXX) -c $(CXXFLAGS) $(LIBS) $< -o $@
%.o: %.c %.h
	$(CC) -c $(CFLAGS) $(LIBS) $< -o $@
%.o: %.f %.h
	$(FC) -c $(FFLAGS) $(LIBS)  $< -o $@ 
%.o: %.f90 
	$(FC) -c $(FFLAGS) $(LIBS) $< -o $@

all: ${objects} ${objectsc}

#find root/a *.c -exec mv {} root/b \;
	
	

clean:
	rm -f ${wildcard *.o}



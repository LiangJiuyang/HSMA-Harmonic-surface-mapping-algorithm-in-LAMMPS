# Harmonic-surface-mapping-algorithm-in-LAMMPS
[HSMA3D](https://aip.scitation.org/doi/10.1063/1.5044438) and [HSMA2D](https://aip.scitation.org/doi/10.1063/5.0003293) (with planar dielectric interfaces) have been implemented into LAMMPS as a k-space module. This module is written via C++ and is paralleled via MPI + OpenMP. We recommend user install 'user-omp' package in Lammps. Fewer MPIs and more OpenMPs are the most efficient choice. We suggest not to use pure MPI. With optimal choice of parameters, the speed of this package is comparable to PPPM (with Intel optimization) in LAMMPS, or even faster than it. This package requires the Intel Parallel Studio or at least the Intel MKL library for solving the least square problem. We offer installation with both `make` and `cmake` (recommended) procedure.

## Installation with Cmake (recommended)
Compiling HSMA using Cmake is a recently-developed project. In this version, one can use a minimal example using the command line version of CMake to build HSMA in LAMMPS is as follows
```
cd lammps                # change to the LAMMPS distribution directory
mkdir build; cd build    # create and use a build directory
cmake -C ../cmake/presets/intel.cmake ../cmake # configuration reading CMake scripts from ../cmake
cmake -D PKG_HSMA=on -D PKG_INTEL=on -D PKG_MOLECULE=on -D PKG_KSPACE=on -D PKG_MANYBODY=on -D INTEL_ARCH=cpu . # include package
cmake --build .          # compilation (or type "make")
```

## Installation with make
For employing HSMA3D after download this full package, the first thing is to include the HSMA package and other appropriate packages in your LAMMPS (cd ./src catalogue):
```
make yes-molecule yes-manybody yes-kspace yes-hsma
```

We recommend to install 'OPENMP' package for the better performance of the evaluation of the LJ potential:
```
make yes-openmp
```

Then, check if following commands is written into the `src/Makefile.package` file : 
```
# Settings for libraries used by specific LAMMPS packages
# this file is auto-edited when those packages are included/excluded
SOURCE = $(wildcard ../HSMA/OFile/*.o)
PKG_INC =   -DLMP_USER_OMP 
PKG_PATH =  $(SOURCE) -lgfortran -ldl
PKG_LIB =   
PKG_CPP_DEPENDS = 
PKG_LINK_DEPENDS = 
PKG_SYSINC =  
PKG_SYSLIB =  
PKG_SYSPATH = 
```

Next, compile the FMM library under ./src/HSMA/OFile catalogue using
```
unzip HSMA1.0.0.zip
make all
```
Note that the default setting is to use `OpenMP` and `gcc`. We also offer command `OMP=OFF` to exclude OpenMP, and `ICC=ON` to replace `gcc` by Intel compiler `icc`, respectively. 

Finally, compile the LAMMPS under ./src catalogue using
```
make intel_cpu_intelmpi
```
which is our recommendation. Other options are also OK like `make omp`, but one should add the dependency `-fopenmp -lmkl_rt` at the end of `PKG_PATH` in `src/Makefile.package` file.

Note that if you want to use HSMA with your own LAMMPS, please just copy /src/HSMA to your lammps/src catalogue, and then refer to the same installation procedure listed in this part.

## User guide
For employing HSMA3D, the only thing needed is to change the k-space solver in your Lammps in-file, just as 
```
kspace_style HSMA3D 1.0e-4 1.2 8 128 55.0 89.0 0 0
```

The first two parameters don't need to modify. Just keep them as `"kspace_style HSMA".`

1.0e-4 : the tolerance of FMM or directly compute.  

1.2 : If the dimensions of simulation box are Lx=Ly=Lz=100, then this parameter means that one find all the images within a sphere with radius 1.2 * sqrt(LxLx+LyLy+LzLz)
.   

8 : means that the number of basis is 8*8=64.  

128 : the number of monitoring points.  

55.0 and 89.0 : the parameters of Fibonacci quadrature. In most cases, 55 and 89 are sufficient so that the relative error of spherical integration is less than 1e−5 (by the rule of thumb). One can set these two numbers as two adjacent Fibonacci numbers, like "89.0 144.0" or "144.0 233.0".

0 (first) : Indicate if one employs FMM (O(N)) to evaluate the potential of near-field. "0" indicates directly computing which has O(N^2) complexity.

0 (second) : Indicate if one employs FMM (O(N)) to evaluate the potential of far-field. "0" indicates directly computing which has O(N^2) complexity.

For employing HSMA2D, the only thing needed is to change the k-space solver in your Lammps in-file, just as 
```
kspace_style HSMA2D 1.0e-3 1.5 0.939 6 40 16 55.0 89.0 1 0
```  
The first two parameters don't need to modify. Just keep them as `"kspace_style HSMA2D".`  

1.0e-3 : the tolerance of FMM or directly compute.  

1.5 : If the dimensions of simulation box are Lx=Ly=Lz=100, then this parameter means that one find all the images within a sphere with radius 1.5 * sqrt(Lx*Lx+Ly*Ly+Lz*Lz)
.  

0.00 : the dielectric mismatch. The range of mismatch is [-1,1].
  
6 : means that the number of basis is 6*6=36. 

40 : the number of Gaussian quadratures. In the current version, this parameter could be selected from {2,3,5,8,10,15,20,30,40,50,60,70,80,90,100,120,160,240,320,480}.   

16 : the parameter "w" in our paper (relative to the 2D dilation quadrature).

55.0 and 89.0 : the parameters of Fibonacci quadrature. In most cases, 55 and 89 are sufficient so that the relative error of spherical integration is less than 1e−5 (by the rule of thumb). One can set these two numbers as two adjacent Fibonacci numbers, like "89.0 144.0" or "144.0 233.0".

1 : Indicate if one employs FMM (O(N)) to evaluate the potential of near-field. "0" indicates directly computing which has O(N^2) complexity.

0 : Indicate if one employs FMM (O(N)) to evaluate the potential of far-field. "0" indicates directly computing which has O(N^2) complexity.

For more details of parameter setting, please refer to our JCP papers which contain the set of parameter within given accuracy. 

Note that the 'pair style' should be set as 'lj/cut' (or lj/cut/omp, we recommend using user-omp package in Lammps) if you want to evaluate LJ potential. Please do not use pair styles which are combined with the near part of a Coulomb solver, such as'lj/cut/coul/long', etc. 

## Examples
Some examples of in-file are provided in the folder `HSMA-Harmonic-surface-mapping-algorithm-in-LAMMPS/HSMA_Example`.
```
salt_1-1.in : 1:1 electrolyte solution (3D)

salt_2-1.in : 2:1 electrolyte solution (3D)

salt_3-1.in : 3:1 electrolyte solution (3D)

2d_salt_2-1.in : 2:1 electrolyte solution (2D with dielectric mismatch)

2d_salt_3-1.in : 3:1 electrolyte solution (2D with dielectric mismatch)
```

To set the number of OpenMP threads per MPI, please type
```
export OMP_NUM_THREADS=40
```
in the command line or dynamiclly set in the code (not recommend). Here `40` should be the number of threads of your machine. Note that this is also set in the beginning of the input script. Please refer to the annotation given in the input files in `HSMA-Harmonic-surface-mapping-algorithm-in-LAMMPS/HSMA_Example`.

To run the input script, an example is cd ./HSMA_Example and then type
```
srun --mpi=pmi2 -n 1 ../src/lmp_intel_cpu_intelmpi -i salt_3-1.in
```
in the command line and hit Enter.

## Citation
In this package, we utilize an efficient implementation of FMM ([FMM3D](https://github.com/flatironinstitute/FMM3D)) which is developed by Greengard's group. 

If you use this package in your work and feel that this package is helpful to you, please cite one (or more) of the following papers in your citation：

1. [Q. Zhao, J. Liang, and Z. Xu, J. Chem. Phys. 149, 084111 (2018).](https://aip.scitation.org/doi/10.1063/1.5044438)
2. [J. Liang, J. Yuan, E. Luijten, and Z. Xu, J. Chem. Phys. 152, 134109 (2020).](https://aip.scitation.org/doi/10.1063/5.0003293)
3. [J. Liang, J. Yuan, and Z. Xu, Comput. Phys. Commun. 276, 108332 (2022).](https://www.sciencedirect.com/science/article/pii/S0010465522000509?via%3Dihub)

This version still need optimization. If you have any questions and suggestions, please send an email to liangjiuyang@sjtu.edu.cn (both Chinese and English are OK).

Good luck to all of you!

## `Introduction of HSMA3D and HSMA2D (with planar dielectric interfaces)`  
Harmonic Surface Mapping Algorithm for 3D periodic systems, first described in paper [Harmonic surface mapping algorithm for fast electrostatic sums](https://aip.scitation.org/doi/10.1063/1.5044438) published by The Journal of Chemical Physics, is an efficient implementation for electrostatic pairwise sums of an infinite number of images accelerated by Fast Multiple method(FMM) and graphics processing units(GPU) (cuda codes are available in another repository). Numerical calculations of the Madelung constant, electrostatic energy of ions in a metallic cavity, and the time performance for large-scale systems show that the HSMA is accurate and fast, and thus is attractive for many applications.

Our recent work, [Harmonic Surface Mapping Algorithm for 2D periodic systems with planar dielectric interfaces](https://aip.scitation.org/doi/10.1063/5.0003293), is published in The Journal of Chemical Physics. We have developed an accurate and efficient method for molecular dynamics simulations of charged particles confined by planar dielectric interfaces. The algorithm combines the image-charge method for near field with the harmonic surface mapping, which converts the contribution of infinite far-field charges into a finite number of charges on an auxiliary spherical surface. We approximate the electrostatic potential of far-field charges via spherical harmonic expansion and determine the coefficients by fitting the Dirichlet-to-Neumann boundary condition, which only requires the potential within the simulation cell. Instead of performing the direct evaluation of spherical harmonic series expansion, we use Green’s second identity to transform the series expansion into a spherical integral, which can be accurately represented by discrete charges on the sphere. Therefore, the fast multipole method can be readily employed to sum over all charges within and on the sphere, achieving truly linear O(N) complexity. Our algorithm can be applied to a broad range of charged complex fluids under dielectric confinement.

Now, we implement two recently developed fast Coulomb solvers, HSMA3D [J. Chem. Phys. 149 (8) (2018) 084111] and HSMA2D [J. Chem. Phys. 152 (13) (2020) 134109], into a new user package [HSMA](https://arxiv.org/abs/2104.05260) for molecular dynamics simulation engine LAMMPS. The HSMA package is designed for efficient and accurate modeling of electrostatic interactions in 3D and 2D periodic systems with dielectric effects at the O(N) cost. The implementation is hybrid MPI and OpenMP parallelized and compatible with existing LAMMPS functionalities. The vectorization technique following AVX512 instructions is adopted for acceleration. To establish the validity of our implementation, we have presented extensive comparisons to the widely used particle-particle particle-mesh (PPPM) algorithm in LAMMPS and other dielectric solvers. With the proper choice of algorithm parameters and parallelization setup, the package enables calculations of electrostatic interactions that outperform the standard PPPM in speed for a wide range of particle numbers.

```
                           Jiuyang Liang
                           Ph.D candidate
                           School of Mathematical Science and Institute of Natural Science
                           Shanghai Jiao Tong University
                           [Homepage in Github](https://github.com/LiangJiuyang/)
                           [Homepage in Researchgate](https://www.researchgate.net/profile/Liang-Jiuyang)
```



<a href="https://info.flagcounter.com/teQT"><img src="https://s01.flagcounter.com/count2/teQT/bg_FFFFFF/txt_000000/border_CCCCCC/columns_5/maxflags_12/viewers_People/labels_0/pageviews_0/flags_0/percent_0/" alt="Flag Counter" border="0"></a>

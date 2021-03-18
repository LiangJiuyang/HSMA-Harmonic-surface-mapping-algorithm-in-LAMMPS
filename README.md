# Harmonic-surface-mapping-algorithm-in-LAMMPS
HSMA3D and HSMA2D (with planar dielectric interfaces) have been implemented into LAMMPS as a k-space module. This module is written via C++ and is paralleled via MPI + OpenMP. We recommend user install 'user-omp' package in Lammps. Fewer MPIs and more OpenMPs are the most efficient choice. We suggest not to use pure MPI.  

For employing HSMA3D, the only thing needed is to change the k-space solver in your Lammps in-file, just as 
```
kspace_style HSMA 1.0e-3 1.3 8 128 55.0 89.0 1 0
```
For employing HSMA2D, the only thing needed is to change the k-space solver in your Lammps in-file, just as 
```
kspace_style HSMA2D 1.0e-3 1.5 0.0 6 40 16 55.0 89.0 1 1
```  

Note that the 'pair style' should be set as 'lj/cut' (or lj/cut/omp, we recommend using user-omp package in Lammps) if you want to evaluate LJ potential. Please do not use pair styles which are combined with the near part of a Coulomb solver, such as'lj/cut/coul/long', etc. 

## `Introduction of HSMA3D and HSMA2D (with planar dielectric interfaces)`  
Harmonic Surface Mapping Algorithm for 3D periodic systems, first described in paper [Harmonic surface mapping algorithm for fast electrostatic sums](https://aip.scitation.org/doi/10.1063/1.5044438) published by The Journal of Chemical Physics, is an efficient implementation for electrostatic pairwise sums of an infinite number of images accelerated by Fast Multiple method(FMM) and graphics processing units(GPU) (cuda codes are available in another repository). Numerical calculations of the Madelung constant, electrostatic energy of ions in a metallic cavity, and the time performance for large-scale systems show that the HSMA is accurate and fast, and thus is attractive for many applications.

Our recent work, [Harmonic Surface Mapping Algorithm for 2D periodic systems with planar dielectric interfaces](https://aip.scitation.org/doi/10.1063/5.0003293), is published in The Journal of Chemical Physics. We have developed an accurate and efficient method for molecular dynamics simulations of charged particles confined by planar dielectric interfaces. The algorithm combines the image-charge method for near field with the harmonic surface mapping, which converts the contribution of infinite far-field charges into a finite number of charges on an auxiliary spherical surface. We approximate the electrostatic potential of far-field charges via spherical harmonic expansion and determine the coefficients by fitting the Dirichlet-to-Neumann boundary condition, which only requires the potential within the simulation cell. Instead of performing the direct evaluation of spherical harmonic series expansion, we use Green’s second identity to transform the series expansion into a spherical integral, which can be accurately represented by discrete charges on the sphere. Therefore, the fast multipole method can be readily employed to sum over all charges within and on the sphere, achieving truly linear O(N) complexity. Our algorithm can be applied to a broad range of charged complex fluids under dielectric confinement.

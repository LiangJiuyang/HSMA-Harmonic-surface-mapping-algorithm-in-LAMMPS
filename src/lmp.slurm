#!/bin/bash

#SBATCH -J lammps 
#SBATCH -p cpu
#SBATCH -N 20   
#SBATCH --ntasks-per-node=40    
#SBATCH --mail-type=end
#SBATCH --mail-user=liangjiuyang@sjtu.edu.cn
#SBATCH -o 800.out
#SBATCH -e %j.err 
#SBATCH --time=04:00:00


module purge
#module load mpich fftw
module load intel-parallel-studio
#module load lammps/20190807-intel-19.0.5-impi 

export I_MPI_PMI_LIBRARY=/usr/lib64/libpmi.so
export I_MPI_FABRICS=shm:ofi  
 
ulimit -s unlimited 
ulimit -l unlimited

srun --mpi=pmi2 ../src/lmp_intel_cpu_intelmpi -i in.spce-bulk-nvt 

mv Time.txt Time_400.txt


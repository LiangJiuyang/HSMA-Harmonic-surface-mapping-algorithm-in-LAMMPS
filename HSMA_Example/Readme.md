#Reference command (for a machine with 40 threads)

#Set the number of OpenMP threads in the environment
export OMP_NUM_THREADS=40

#run with 'srun'
srun --mpi=pmi2 -n 1 ../src/lmp_intel_cpu_intelmpi -i salt_3-1.in
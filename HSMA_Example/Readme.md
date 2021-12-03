#Reference command (for a machine with 40 threads)

#Set the number of OpenMP threads in the environment
```
export OMP_NUM_THREADS=40
```

#run with `srun` with `make` for installing
```
srun --mpi=pmi2 -n 1 ../src/lmp_intel_cpu_intelmpi -i salt_3-1.in
```

#Use cmake for installing
```
srun --mpi=pmi2 -n 1 ../build/lmp -i salt_3-1.in
```

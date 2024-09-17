# nBody
This project implements a brute-force N-body gravitational simulation in C. The simulation calculates gravitational interactions between a large number of bodies over several iterations. It demonstrates performance improvements achieved using parallel processing techniques, specifically MPI (Message Passing Interface) and OpenMP.

Project Overview:

Serial Implementation: Calculates gravitational forces between all pairs of bodies using a nested loop. Serves as a baseline for performance comparison.
MPI Parallelization: Parallelizes the simulation using MPI to distribute the workload across multiple ranks, significantly reducing runtime.
OpenMP Optimization: Uses OpenMP to further parallelize the simulation on shared-memory systems for multi-core CPUs.
Performance Summary:

Serial Implementation: Average runtime of 10 seconds per iteration with 93 million interactions per second.
MPI Parallel Implementation with 4 ranks: Average runtime of 3.5 seconds per iteration with 220 million interactions per second.
Requirements:

- C Compiler (e.g., gcc)
- MPI Library (e.g., Open MPI)
- OpenMP support (included in most modern C compilers)
Usage:

Serial Implementation: Compile with "gcc -o nbody_serial nbody_serial.c -lm" and run with "./nbody_serial"
MPI Implementation: Compile with "mpicc -o nbody_mpi nbody_mpi.c -lm" and run with "mpirun -np 4 ./nbody_mpi"
OpenMP Implementation: Compile with "gcc -o nbody_openmp nbody_openmp.c -lm -fopenmp" and run with "export OMP_NUM_THREADS=4; ./nbody_openmp"
This project demonstrates the use of parallel processing to significantly reduce computation time in an N-body gravitational simulation. It highlights the application of both MPI for distributed-memory parallelism and OpenMP for shared-memory parallelism.

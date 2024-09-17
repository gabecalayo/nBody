N-Body Gravitational Simulation
This project implements a brute-force N-body gravitational simulation in C. The simulation calculates the gravitational interactions between a large number of bodies over several iterations. It is designed to demonstrate the performance improvements achieved by leveraging parallel processing techniques, specifically using MPI (Message Passing Interface) and OpenMP.

Project Overview
Serial Implementation: A brute-force approach that calculates the gravitational forces between all pairs of bodies. This implementation serves as a baseline for performance comparison.
MPI Parallelization: The simulation is parallelized using MPI, enabling distributed processing across multiple ranks. This approach significantly reduces runtime by dividing the workload among processes.
OpenMP Optimization: Additionally, the simulation can be parallelized using OpenMP for shared-memory systems, offering further performance enhancements on multi-core CPUs.
Features
Serial Simulation: Computes gravitational interactions for a given number of bodies using a nested loop, resulting in O(n²) complexity.
MPI Parallelization: Distributes the computational workload across multiple processes, reducing runtime and increasing the number of interactions per second.
OpenMP Parallelization: Optimizes the simulation for multi-core processors, improving performance on shared-memory systems.
Performance Metrics: Measures and reports runtime and interactions per second for different implementations, allowing for detailed performance analysis.
Performance Summary
Serial Implementation:

Average runtime: 10 seconds per iteration
Interaction rate: 93 million interactions/second
MPI Parallel Implementation (with 4 ranks):

Average runtime: 3.5 seconds per iteration
Interaction rate: 220 million interactions/second
Requirements
C Compiler (e.g., gcc)
MPI Library (e.g., Open MPI)
OpenMP (Supported by most modern C compilers)
Usage
Compilation
Serial Implementation:

bash
Copy code
gcc -o nbody_serial nbody_serial.c -lm
MPI Implementation:

bash
Copy code
mpicc -o nbody_mpi nbody_mpi.c -lm
OpenMP Implementation:

bash
Copy code
gcc -o nbody_openmp nbody_openmp.c -lm -fopenmp
Running the Simulation
Serial Implementation:

bash
Copy code
./nbody_serial
MPI Implementation (with 4 processes):

bash
Copy code
mpirun -np 4 ./nbody_mpi
OpenMP Implementation:

bash
Copy code
export OMP_NUM_THREADS=4  # Set the number of threads
./nbody_openmp
Example Output
The program will output the runtime for each iteration and the average number of interactions per second:

yaml
Copy code
Iteration 1: 10.000 seconds
Iteration 2: 3.500 seconds
...
Average time per iteration: 3.500 seconds
Performance: 0.220 Billion Interactions / second
Project Structure
nbody_serial.c: The serial implementation of the N-body simulation.
nbody_mpi.c: The MPI-parallelized implementation.
nbody_openmp.c: The OpenMP-parallelized implementation.
timer.h: Header file containing functions for timing execution.
Key Learnings
Parallel Computing: Demonstrates the use of MPI for distributed memory parallelism and OpenMP for shared memory parallelism.
Performance Optimization: Shows how parallel processing can significantly reduce computation time and improve efficiency in computationally intensive tasks.
License
This project is open-source and free to use under the MIT License.

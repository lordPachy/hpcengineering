# OpenMP

We can compile with
`g++ -fopenmp <filename>`

and we can decide the number of threads at execution time with the variable
`OMP_NUM_THREADS=<number of threads>`

## Critical sectons
We can use `#pragma omp critical`.

# MPI

The OpenMPI v5.0 version reference is the indicated one. There is also\

https://docs.open-mpi.org/en/v5.0.9rc2/man-openmpi/man3/index.html\

mpirun and mpiexec are equivalent in OpenMPI.\

Note that when we run things interactively from the cluster shell, only one node is considered. If we want to use more than one, we should pass through PBS.

## Communicators

We can use a non-global communicator to diminish the scope of our communication.
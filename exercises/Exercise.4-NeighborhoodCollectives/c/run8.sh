#!/bin/bash
mpirun -np 8 ./stencil_mpi_carttopo_neighcolls 512 10 10000 4 2

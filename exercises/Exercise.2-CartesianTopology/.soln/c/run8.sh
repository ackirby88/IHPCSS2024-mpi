#!/bin/bash
mpirun -np 8 ./mpi_stencil_cart 512 10 10000 4 2

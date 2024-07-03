#include "mpi.h"
#include <stdio.h>
#define NELEM 25

int main(int argc, char *argv[])  {
    int numtasks, rank, source=0, dest, tag=1, i;

    typedef struct {
        float x, y, z;
        float velocity;
        int  n, type;
    } Particle;

    Particle p[NELEM], particles[NELEM];
    MPI_Datatype particletype, oldtypes[2];   // required variables
    int blockcounts[2];

    // MPI_Aint type used to be consistent with syntax of
    // MPI_Type_extent routine
    MPI_Aint offsets[2], extent;
    MPI_Status stat;

    MPI_Init(&argc,&argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &numtasks);

    /* ===================================================================== */
    /* Step 1. Create an MPI Struct Type
     *    Summary:
     *      Make a new struct derived datatype.
     *
     *    Function Call:
     *       int MPI_Type_struct(int count,
     *                           const int *array_of_blocklengths,
     *                           const MPI_Aint *array_of_displacements,
     *                           const MPI_Datatype *array_of_types,
     *                           MPI_Datatype *newtype);
     *
     *   Input Parameters:
     *       count
     *          number of blocks (integer) -- also number of entries in arrays array_of_types , array_of_displacements and array_of_blocklengths
     *       array_of_blocklengths
     *          number of elements in each block (array)
     *       array_of_displacements
     *          byte displacement of each block (array)
     *       array_of_types
     *          type of elements in each block (array of handles to datatype objects)
     *
     *   Output Parameters:
     *      newtype
     *          new datatype (handle)
     */
    // TODO: setup description of the 4 MPI_FLOAT fields x, y, z, velocity
    offsets[0] = //TODO;
    oldtypes[0] = //TODO;
    blockcounts[0] = //TODO;

    // TODO: setup description of the 2 MPI_INT fields n, type.
    //      We first figure offset by getting size of MPI_FLOAT
    MPI_Type_extent(MPI_FLOAT, &extent);
    offsets[1] = //TODO; HINT: how many 'extent's do we need?
    oldtypes[1] = //TODO;
    blockcounts[1] = //TODO;
    
    // TODO: create the struct data type

    // TODO: commit the new derived datatype 

    /* ===================================================================== */

    // task 0 initializes the particle array and then sends it to each task
    if (rank == 0) {
        for (i=0; i<NELEM; i++) {
            particles[i].x = i * 1.0;
            particles[i].y = i * -1.0;
            particles[i].z = i * 1.0; 
            particles[i].velocity = 0.25;
            particles[i].n = i;
            particles[i].type = i % 2; 
        }
        for (i=0; i<numtasks; i++) {
            MPI_Send(particles, NELEM, particletype, i, tag, MPI_COMM_WORLD);
        }
    }

    // all tasks receive particletype data
    MPI_Recv(p, NELEM, particletype, source, tag, MPI_COMM_WORLD, &stat);

    printf("rank= %d   %3.2f %3.2f %3.2f %3.2f %d %d\n", rank,p[3].x,
           p[3].y,p[3].z,p[3].velocity,p[3].n,p[3].type);

    // free datatype when done using it
    MPI_Type_free(&particletype);
    MPI_Finalize();
    return 0;
}

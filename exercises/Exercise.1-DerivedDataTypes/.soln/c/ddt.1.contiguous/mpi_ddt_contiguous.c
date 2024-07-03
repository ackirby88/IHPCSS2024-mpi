#include "mpi.h"
#include <stdio.h>
#define SIZE 4

int main(int argc, char *argv[])  {
    int numtasks, rank, source=0, dest, tag=1, i;
    float a[SIZE][SIZE] =
        { 1.0, 2.0, 3.0, 4.0,
          5.0, 6.0, 7.0, 8.0,
          9.0,10.0,11.0,12.0,
         13.0,14.0,15.0,16.0};
    float b[SIZE];

    MPI_Status stat;
    MPI_Datatype rowtype;   // required variable

    MPI_Init(&argc,&argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &numtasks);

    /* ===================================================================== */
    /* Step 1. Create an MPI Contiguous Type
     *    Summary:
     *      Make a new contiguous derived datatype.
     *
     *    Function Call:
     *      int MPI_Type_contiguous(int count,
     *                              MPI_Datatype oldtype,
     *                              MPI_Datatype *newtype);
     *      
     *   Input Parameters:
     *       count
     *           replication count (non-negative integer)
     *       oldtype
     *           old datatype (handle)
     *
     *   Output Parameters
     *       newtype
     *           new datatype (handle)
     */
    // TODO: create the contiguous data type
    MPI_Type_contiguous(SIZE, MPI_FLOAT, &rowtype);
    // TODO: commit the new derived datatype 
    MPI_Type_commit(&rowtype);
    /* ===================================================================== */   

    if (numtasks == SIZE) {
        // task 0 sends one element of rowtype to all tasks
        if (rank == 0) {
            /* =================================================================== */
            /* Step 2. Send continguous data type using MPI_Send.
              *    Summary:
              *      Call MPI_Send 
              *
              *    Function Call:
              *      int MPI_Send(const void *buf,
              *                   int count,
              *                   MPI_Datatype datatype,
              *                   int dest,
              *                   int tag,
              *                   MPI_Comm comm);
              *      
              *   Input Parameters:
              *     buf
              *         initial address of send buffer (choice)
              *     count
              *         number of elements in send buffer (non-negative integer)
              *     datatype
              *         datatype of each send buffer element (handle)
              *     dest
              *         rank of destination (integer)
              *     tag
              *         message tag (integer)
              *     comm
              *         communicator (handle)
              */
            // TODO: send each ROW i of the array 'a' using the derived data type.
            for (i=0; i<numtasks; i++) {
                MPI_Send(&a[i][0], 1, rowtype, i, tag, MPI_COMM_WORLD);
            }
            /* =================================================================== */            
        }

        // all tasks receive rowtype data from task 0
        MPI_Recv(b, SIZE, MPI_FLOAT, source, tag, MPI_COMM_WORLD, &stat);
        printf("rank= %d  b= %3.1f %3.1f %3.1f %3.1f\n",
               rank,b[0],b[1],b[2],b[3]);
    } else {
        printf("Must specify %d processors. Terminating.\n",SIZE);
    }

    // free datatype when done using it
    MPI_Type_free(&rowtype);
    MPI_Finalize();
    return 0;
}

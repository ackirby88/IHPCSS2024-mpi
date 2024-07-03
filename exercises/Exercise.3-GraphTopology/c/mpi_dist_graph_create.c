#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>

//number of cpu requested
#define ncpu 4

/* Demonstration of MPI_Dist_graph_create 
 *  with rank reordering.  
 */
int main(int argc, char *argv[]){
    int numtasks, rank, rankG;

    MPI_Init(&argc,&argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &numtasks);
    
    if (numtasks != ncpu) {
        printf("Must specify %d processors. Terminating.\n",ncpu);
	    return 0;
    }

    //Setting up of the graph_comm
    int source[2] = {rank,rank};
    int degree[2];
    int dest[ncpu];
    int weight[ncpu] = {1, 1, 1, 1};
    int recv[ncpu] = {-1, -1, -1, -1};
    int send = rank;

    //Hardcoding the edges: 0<->2, 1<->3
    if (rank == 0) {
        dest[0] = 2;
        dest[1] = 2;
        degree[0] = 1;
        degree[1] = 1;
    } else
    if (rank == 1) {
        dest[0] = 3;
        dest[1] = 3;
        degree[0] = 1;
        degree[1] = 1;
    } else
    if (rank == 2) {
        dest[0] = 0;
        dest[1] = 0;
        degree[0] = 1;
        degree[1] = 1;
    } else
    if (rank == 3) {
        dest[0] = 1;
        dest[1] = 1;
        degree[0] = 1;
        degree[1] = 1;
    }

    /* ===================================================================== */
    /* Step 1. Create an MPI Dist graph.
     *    Summary:
     *      Make a new contiguous derived datatype.
     *
     * int MPI_Dist_graph_create(MPI_Comm comm_old, int n, const int sources[],
     *                           const int degrees[], const int destinations[],
     *                           const int weights[],
     *                           MPI_Info info, int reorder, MPI_Comm *comm_dist_graph);
     *      
     *   Input Parameters:
     *       comm_old
     *          input communicator (handle)
     *       n
     *          number of source nodes for which this process specifies edges (non-negative integer)
     *       sources
     *          array containing the n source nodes for which this process specifies edges (array of non-negative integers)
     *       degrees
     *          array specifying the number of destinations for each source node in the source node array (array of non-negative integers)
     *       destinations
     *          destination nodes for the source nodes in the source node array (array of non-negative integers)
     *       weights
     *          weights for source to destination edges (array of non-negative integers or MPI_UNWEIGHTED)
     *       info
     *          hints on optimization and interpretation of weights (handle)
     *       reorder
     *          the process may be reordered (true) or not (false) (logical)
     *  Output Parameters
     *       comm_dist_graph
     *          communicator with distributed graph topology added (handle)
     */
    MPI_Comm graph_comm;
    // TODO: create the distributed graph communicator (graph_comm)
    
    // TODO: get the new rank ID (rankG) in the new communicator
    
    /* ===================================================================== */

    //print the result
    char name[MPI_MAX_PROCESSOR_NAME];
    int resultlen;
    MPI_Get_processor_name(name, &resultlen);

    printf("Rank: %i on %s, Rank reordered: %i \n", rank, name, rankG);

    //-----------------------------------------------
    // Print the communication graph_comm, just to be sure
    MPI_Barrier(MPI_COMM_WORLD);

    /* ===================================================================== */
    /* Step 2. Get neighbors in new distributed graph.
     *    Summary:
     *      Make a new contiguous derived datatype.
     *
     * ---------------------------------------------------------------------
     * int MPI_Dist_graph_neighbors_count(MPI_Comm comm, int *indegree, 
     *                                    int *outdegree, int *weighted);
     *      
     *   Input Parameters:
     *       comm
     *           communicator with distributed graph topology (handle)
     *
     *   Output Parameters:
     *       indegree
     *           number of edges into this process (non-negative integer)
     *       outdegree
     *           number of edges out of this process (non-negative integer)
     *       weighted
     *           false if MPI_UNWEIGHTED was supplied during creation, 
     *           true otherwise (logical)
     * ---------------------------------------------------------------------
     * int MPI_Dist_graph_neighbors(MPI_Comm comm,
     *                              int maxindegree, int sources[], int sourceweights[],
     *                              int maxoutdegree, int destinations[], int destweights[]);
     *
     *   Input Parameters:
     *       comm
     *           communicator with distributed graph topology (handle)
     *       maxindegree
     *           size of sources and sourceweights arrays (non-negative integer)
     *       maxoutdegree
     *           size of destinations and destweights arrays (non-negative integer)
     *
     *   Output Parameters:
     *       sources
     *           processes for which the calling process is a destination (array of non-negative integers)
     *       sourceweights
     *           weights of the edges into the calling process (array of non-negative integers)
     *       destinations
     *           processes for which the calling process is a source (array of non-negative integers)
     *       destweights
     *           weights of the edges out of the calling process (array of non-negative integers)
     * ---------------------------------------------------------------------
     */
    int inD, outD, wei;
    // TODO: get the number of graph neighbors (store in variables: inD, outD, wei)
    
    printf("IN-Degree:%d OUT-Degree:%d weight:%d\n",inD,outD,wei);

    // allocating the source and destination arrays based on counts
    int* Sour = (int*) malloc(sizeof(int)*inD);
    int* SourW = (int*) malloc(sizeof(int)*inD);
    int* Dest = (int*) malloc(sizeof(int)*outD);
    int* DestW = (int*) malloc(sizeof(int)*outD);

    // TODO: get the graph neighbors
    
    /* ===================================================================== */
    
    printf("IN-Edges: ");
    for(int i=0; i<inD; i++) printf("%d (wght:%d), ",Sour[i],SourW[i]);
    
    printf("\n OUT-Edges: ");
    for(int i=0; i<outD; i++) printf("%d (wght:%d), ",Dest[i],DestW[i]);
    printf("\n");
    
    MPI_Finalize();
    return 0;
}

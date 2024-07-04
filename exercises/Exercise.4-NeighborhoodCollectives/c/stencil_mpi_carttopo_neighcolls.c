#include "stencil_par.h"

int main(int argc, char **argv) {

  MPI_Init(&argc, &argv);
  int r,p;
  MPI_Comm comm = MPI_COMM_WORLD;
  MPI_Comm_rank(comm, &r);
  MPI_Comm_size(comm, &p);
  int n, energy, niters;

  if (r==0) {
     // argument checking
      if(argc < 4) {
          if(!r) printf("usage: stencil_mpi <n> <energy> <niters>\n");
          MPI_Finalize();
          exit(1);
      }

      n = atoi(argv[1]); // nxn grid
      energy = atoi(argv[2]); // energy to be injected per iteration
      niters = atoi(argv[3]); // number of iterations

      // distribute arguments
      int args[3] = {n, energy, niters};
      MPI_Bcast(args, 3, MPI_INT, 0, comm);
  }
  else {
      int args[3];
      MPI_Bcast(args, 3, MPI_INT, 0, comm);
      n=args[0]; energy=args[1]; niters=args[2];
  }

  int pdims[2]={0,0};
  // compute good (rectangular) domain decomposition
  MPI_Dims_create(p, 2, pdims);
  int px = pdims[0];
  int py = pdims[1];

  // create Cartesian topology
  int periods[2] = {0,0};
  MPI_Comm topocomm;
  MPI_Cart_create(comm, 2, pdims, periods, 0, &topocomm);

  // get my local x,y coordinates
  int coords[2];
  MPI_Cart_coords(topocomm, r, 2, coords);
  int rx = coords[0];
  int ry = coords[1];

  int source, north, south, east, west;
  MPI_Cart_shift(topocomm, 0, 1, &west, &east);
  MPI_Cart_shift(topocomm, 1, 1, &north, &south);

  // decompose the domain
  int bx = n/px; // block size in x
  int by = n/py; // block size in y
  int offx = rx*bx; // offset in x
  int offy = ry*by; // offset in y

  //  printf("%i (%i,%i) - w: %i, e: %i, n: %i, s: %i\n", r, ry,rx,west,east,north,south);

  double *aold = (double*)calloc(1,(bx+2)*(by+2)*sizeof(double)); // 1-wide halo zones!
  double *anew = (double*)calloc(1,(bx+2)*(by+2)*sizeof(double)); // 1-wide halo zones!
  double *tmp;

  #define nsources 3
  int sources[nsources][2] = {{n/2,n/2}, {n/3,n/3}, {n*4/5,n*8/9}};
  int locnsources=0; // number of sources in my area
  int locsources[nsources][2]; // sources local to my rank
  for (int i=0; i<nsources; ++i) { // determine which sources are in my patch
    int locx = sources[i][0] - offx;
    int locy = sources[i][1] - offy;
    if(locx >= 0 && locx < bx && locy >= 0 && locy < by) {
      locsources[locnsources][0] = locx+1; // offset by halo zone
      locsources[locnsources][1] = locy+1; // offset by halo zone
      locnsources++;
    }
  }
  //printf("%i: locnsources: %i\n", r, locnsources);

  double heat; // total heat in system
  for(int iter=0; iter<niters; ++iter) {
    // refresh heat sources
    for(int i=0; i<locnsources; ++i) {
      aold[ind(locsources[i][0],locsources[i][1])] += energy; // heat source
    }

    // exchange data with neighbors
    double *sbuf = (double*)calloc(1,2*bx*sizeof(double)+2*by*sizeof(double)); // send buffer (west, east, north, south)
    double *rbuf = (double*)calloc(1,2*bx*sizeof(double)+2*by*sizeof(double)); // receive buffer (w, e, n, s)
    for(int i=0; i<by; ++i) sbuf[i] = aold[ind(1,i+1)]; // pack west
    for(int i=0; i<by; ++i) sbuf[by+i] = aold[ind(bx,i+1)]; // pack east
    for(int i=0; i<bx; ++i) sbuf[2*by+i] = aold[ind(i+1,1)]; // pack north
    for(int i=0; i<bx; ++i) sbuf[2*by+bx+i] = aold[ind(i+1,by)]; // pack south

    /* ===================================================================== */
    /* Step 1. Create the neighborhood collective call.
     *    Summary:
     *      Perform a nonblocking neighborhood collective call using variable
     *      data sizes.
     *
     *    Function Call:
     *      int MPI_Ineighbor_alltoallv(const void *sendbuf, const int sendcounts[], const int sdispls[], MPI_Datatype sendtype,
     *                                  void *recvbuf, const int recvcounts[], const int rdispls[], MPI_Datatype recvtype,
     *                                  MPI_Comm comm, MPI_Request *request);
     *    Input Parameters:
     *        sendbuf
     *            starting address of the send buffer (choice)
     *        sendcounts
     *            non-negative integer array (of length outdegree) specifying the number of elements to send to each neighbor
     *        sdispls
     *            integer array (of length outdegree). Entry j specifies the displacement (relative to sendbuf) from which to send the outgoing data to neighbor j
     *        sendtype
     *            data type of send buffer elements (handle)
     *        recvcounts
     *            non-negative integer array (of length indegree) specifying the number of elements that are received from each neighbor
     *        rdispls
     *            integer array (of length indegree). Entry i specifies the displacement (relative to recvbuf) at which to place the incoming data from neighbor i.
     *        recvtype
     *            data type of receive buffer elements (handle)
     *        comm
     *            communicator with topology structure (handle)
     *
     *    Output Parameters:
     *        recvbuf
     *            starting address of the receive buffer (choice)
     *        request
     *            communication request (handle)
     */
    // TODO: fillout 'counts' array: how much elements of data are we sending to each neighbor (hint: look at the packing of sbuf above)
    int counts[4] = {//TODO, //TODO, //TODO, //TODO};
    // TODO: fillout 'displs' array: What are the starting indexes of each "chunk" (hint: look at the packing of sbuf above)
    int displs[4] = {//TODO, //TODO, //TODO, //TODO};

    MPI_Request req;
    MPI_Status status;
    // TODO: perform nonblocking neighborhood collective call (note: counts serves for both sendcounts and recvcounts; displs serves for both sdispls and rdispls)
    
    MPI_Wait(&req, &status);
    /* ===================================================================== */

    for(int i=0; i<by; ++i) aold[ind(0,i+1)] = rbuf[i]; // unpack loop
    for(int i=0; i<by; ++i) aold[ind(bx+1,i+1)] = rbuf[by+i]; // unpack loop
    for(int i=0; i<bx; ++i) aold[ind(i+1,0)] = rbuf[2*by+i]; // unpack loop
    for(int i=0; i<bx; ++i) aold[ind(i+1,by+1)] = rbuf[2*by+bx+i]; // unpack loop

    heat = 0.0;
    for(int j=1; j<by+1; ++j) {
      for(int i=1; i<bx+1; ++i) {
        anew[ind(i,j)] = aold[ind(i,j)]/2.0 + (aold[ind(i-1,j)] + aold[ind(i+1,j)] + aold[ind(i,j-1)] + aold[ind(i,j+1)])/4.0/2.0;
        heat += anew[ind(i,j)];
      }
    }

    tmp=anew; anew=aold; aold=tmp; // swap arrays
    if(iter == niters-1) printarr_par(iter, anew, n, px, py, rx, ry, bx, by, offx, offy, comm);
  }
  double rheat;
  MPI_Allreduce(&heat, &rheat, 1, MPI_DOUBLE, MPI_SUM, comm);
  if(!r) printf("[%i] last heat: %f\n", r, rheat);

  MPI_Finalize();
}

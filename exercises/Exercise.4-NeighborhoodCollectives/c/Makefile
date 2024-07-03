MPICC=mpicc
CFLAGS=-std=c99
LIBS=-lm
OBJS=printarr_par.o
BINS=stencil_mpi_carttopo_neighcolls

all: $(BINS)

%.o: %.c Makefile
	$(MPICC) $(CFLAGS) $< -c -o $@

stencil_mpi_carttopo_neighcolls: stencil_mpi_carttopo_neighcolls.c Makefile $(OBJS) 
	$(MPICC) $(CFLAGS) -o $@ $< $(OBJS) $(LIBS)

clean:
	rm -f $(BINS) $(OBJS)
	rm -f output*.bmp

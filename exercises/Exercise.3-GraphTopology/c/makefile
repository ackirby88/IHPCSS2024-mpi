CXX = mpicxx
CFLAGS = -g
EXEC = mpi_dist_graph_create

all:
	$(CXX) $(CFLAGS) mpi_dist_graph_create.c -o $(EXEC)
clean:
	@echo "*** Clean : removing object files"
	@( $(RM) *.o *.mod *~)
	@( $(RM) $(EXEC))

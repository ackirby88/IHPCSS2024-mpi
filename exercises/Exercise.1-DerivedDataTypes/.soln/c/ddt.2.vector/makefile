CXX = mpicxx
CFLAGS = -g
EXEC = mpi_ddt_vector

all:
	$(CXX) $(CFLAGS) mpi_ddt_vector.c -o $(EXEC)
clean:
	@echo "*** Clean : removing object files"
	@( $(RM) *.o *.mod *~)
	@( $(RM) $(EXEC))

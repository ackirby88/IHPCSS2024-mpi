CXX = mpicxx
CFLAGS = -g
EXEC = mpi_ddt_contiguous

all:
	$(CXX) $(CFLAGS) mpi_ddt_contiguous.c -o $(EXEC)
clean:
	@echo "*** Clean : removing object files"
	@( $(RM) *.o *.mod *~)
	@( $(RM) $(EXEC))

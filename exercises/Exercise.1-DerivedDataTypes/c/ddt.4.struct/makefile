CXX = mpicxx
CFLAGS = -g
EXEC = mpi_ddt_struct

all:
	$(CXX) $(CFLAGS) mpi_ddt_struct.c -o $(EXEC)
clean:
	@echo "*** Clean : removing object files"
	@( $(RM) *.o *.mod *~)
	@( $(RM) $(EXEC))

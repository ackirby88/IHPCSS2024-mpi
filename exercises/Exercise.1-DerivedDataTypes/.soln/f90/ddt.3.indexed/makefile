FC = mpif90
FFLAGS = -g
EXEC = mpi_ddt_indexed

all:
	$(FC) $(FFLAGS) mpi_ddt_indexed.f90 -o $(EXEC)
clean:
	@echo "*** Clean : removing object files"
	@( $(RM) *.o *.mod *~)
	@( $(RM) $(EXEC))

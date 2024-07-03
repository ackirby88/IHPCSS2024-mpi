FC = mpif90
FFLAGS = -g
EXEC = mpi_ddt_struct

all:
	$(FC) $(FFLAGS) mpi_ddt_struct.f90 -o $(EXEC)
clean:
	@echo "*** Clean : removing object files"
	@( $(RM) *.o *.mod *~)
	@( $(RM) $(EXEC))

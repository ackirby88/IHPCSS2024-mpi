FC = mpif90
CFLAGS = -g
EXEC = mpi_ddt_contiguous

all:
	$(FC) $(FFLAGS) mpi_ddt_contiguous.f90 -o $(EXEC)
clean:
	@echo "*** Clean : removing object files"
	@( $(RM) *.o *.mod *~)
	@( $(RM) $(EXEC))

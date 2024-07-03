FC = mpif90
fFLAGS = -g
EXEC = mpi_ddt_vector

all:
	$(FC) $(FFLAGS) mpi_ddt_vector.f90 -o $(EXEC)
clean:
	@echo "*** Clean : removing object files"
	@( $(RM) *.o *.mod *~)
	@( $(RM) $(EXEC))

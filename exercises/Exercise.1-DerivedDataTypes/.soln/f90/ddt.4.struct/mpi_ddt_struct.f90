program struct
    include 'mpif.h'

    integer NELEM
    parameter(NELEM=25)
    integer numtasks, rank, source, dest, tag, i,  ierr
    integer stat(MPI_STATUS_SIZE)

    type Particle
      sequence
        real*4 x, y, z, velocity
        integer n, type
    end type Particle

    type (Particle) p(NELEM), particles(NELEM)
    integer particletype, oldtypes(0:1)   ! required variables
    integer blockcounts(0:1), offsets(0:1), extent
    tag = 1

    call MPI_INIT(ierr)
    call MPI_COMM_RANK(MPI_COMM_WORLD, rank, ierr)
    call MPI_COMM_SIZE(MPI_COMM_WORLD, numtasks, ierr)

    ! =====================================================================
    ! Step 1. Create an MPI Struct Type
    !    Summary:
    !      Make a new struct derived datatype.
    !
    !    Function Call:
    !       MPI_TYPE_STRUCT(INTEGER COUNT, 
    !                       INTEGER(*) ARRAY_OF_BLOCKLENGTHS,
    !                       INTEGER(*) ARRAY_OF_DISPLACEMENTS,
    !                       INTEGER(*) ARRAY_OF_TYPES,
    !                       INTEGER NEWTYPE,
    !                       INTEGER IERROR)
    !
    !   Input Parameters:
    !       COUNT
    !          number of blocks (integer) -- also number of entries in arrays array_of_types , array_of_displacements and array_of_blocklengths
    !       ARRAY_OF_BLOCKLENGTHS
    !          number of elements in each block (array)
    !       ARRAY_OF_DISPLACEMENTS
    !          byte displacement of each block (array)
    !       ARRAY_OF_TYPES
    !          type of elements in each block (array of handles to datatype objects)
    !
    !   Output Parameters:
    !      NEWTYPE
    !          new datatype (handle)
    !      IERROR
    !          Fortran only: Error status (integer).
    !
    ! TODO: setup description of the 4 MPI_REAL fields x, y, z, velocity
    offsets(0) = 0 ! TODO
    oldtypes(0) = MPI_REAL ! TODO
    blockcounts(0) = 4 ! TODO

    ! TODO: setup description of the 2 MPI_INTEGER fields n, type 
    !       We first figure offset by getting size of MPI_REAL
    call MPI_TYPE_EXTENT(MPI_REAL, extent, ierr)
    offsets(1) = 4 * extent ! TODO; HINT: how many 'extent's do we need?
    oldtypes(1) = MPI_INTEGER  ! TODO 
    blockcounts(1) = 2  ! TODO

    ! TODO: create the struct data type
    call MPI_TYPE_STRUCT(2, blockcounts, offsets, oldtypes, particletype, ierr)
    ! TODO: commit the new derived datatype 
    call MPI_TYPE_COMMIT(particletype, ierr)
    ! =====================================================================

    ! task 0 initializes the particle array and then sends it to each task
    if (rank .eq. 0) then
        do i=1, NELEM
            particles(i) = Particle ( 1.0*i, -1.0*i, 1.0*i, 0.25, i, mod(i,2) )
        end do

        do i=0, numtasks-1
            call MPI_SEND(particles, NELEM, particletype, i, tag, &
                          MPI_COMM_WORLD, ierr)
        end do
    endif

    ! all tasks receive particletype data
    source = 0
    call MPI_RECV(p, NELEM, particletype, source, tag, &
                  MPI_COMM_WORLD, stat, ierr)

    print *, 'rank= ',rank,' p(3)= ',p(3)

    ! free datatype when done using it
    call MPI_TYPE_FREE(particletype, ierr)
    call MPI_FINALIZE(ierr)
end program

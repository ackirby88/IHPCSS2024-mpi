program indexed
    include 'mpif.h'

    integer NELEMENTS
    parameter(NELEMENTS=6)
    integer numtasks, rank, source, dest, tag, i,  ierr
    integer blocklengths(0:1), displacements(0:1)
    real*4 a(0:15), b(0:NELEMENTS-1)
    integer stat(MPI_STATUS_SIZE)
    integer indextype   ! required variable
    tag = 1

    data a  /1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, &
            9.0, 10.0, 11.0, 12.0, 13.0, 14.0, 15.0, 16.0 /

    call MPI_INIT(ierr)
    call MPI_COMM_RANK(MPI_COMM_WORLD, rank, ierr)
    call MPI_COMM_SIZE(MPI_COMM_WORLD, numtasks, ierr)

    ! =====================================================================
    ! Step 1. Create an MPI Indexed Type
    !    Summary:
    !      Make a new indexed derived datatype.
    !
    !    Function Call:
    !      MPI_TYPE_INDEXED(INTEGER COUNT,
    !                       INTEGER(*) ARRAY_OF_BLOCKLENGTHS,
    !                       INTEGER(*) ARRAY_OF_DISPLACEMENTS,
    !                       INTEGER OLDTYPE,
    !                       INTEGER NEWTYPE,
    !                       INTEGER IERROR)
    !      
    !   Input Parameters:
    !      COUNT
    !          number of blocks -- also number of entries in array_of_displacements and array_of_blocklengths
    !      ARRAY_OF_BLOCKLENGTHS
    !          number of elements in each block (array of nonnegative integers)
    !      ARRAY_OF_DISPLACEMENTS
    !          displacement of each block in multiples of oldtype (array of integers)
    !      OLDTYPE
    !          old datatype (handle)
    !
    !   Output Parameters
    !      NEWTYPE
    !          new datatype (handle)
    !      IERROR
    !          Fortran only: Error status (integer).
    !
    ! NOTE: We want the the resulting values of b[] to be {6.0 7.0 8.0 9.0 13.0 14.0}.
    ! TODO: fill in the values for blocklengths
    blocklengths(0) = 4 ! TODO
    blocklengths(1) = 2 ! TODO

    ! TODO: fill in the values for the displacements
    displacements(0) = 5  ! TODO
    displacements(1) = 12 ! TODO
    
    ! TODO: create the indexed data type
    call MPI_TYPE_INDEXED(2, blocklengths, displacements, MPI_REAL, indextype, ierr)
    ! TODO: commit the new derived datatype 
    call MPI_TYPE_COMMIT(indextype, ierr)
    ! =====================================================================

    if (rank .eq. 0) then
        ! task 0 sends one element of indextype to all tasks
        do i=0, numtasks-1
            call MPI_SEND(a, 1, indextype, i, tag, MPI_COMM_WORLD, ierr)
        end do
    endif

    ! all tasks receive indextype data from task 0
    source = 0
    call MPI_RECV(b, NELEMENTS, MPI_REAL, source, tag, MPI_COMM_WORLD, &
                    stat, ierr)
    print *, 'rank= ',rank,' b= ',b

    ! free datatype when done using it
    call MPI_TYPE_FREE(indextype, ierr)
    call MPI_FINALIZE(ierr)

end program
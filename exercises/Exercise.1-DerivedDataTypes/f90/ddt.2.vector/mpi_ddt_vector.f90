program vector
    include 'mpif.h'

    integer SIZE
    parameter(SIZE=4)
    integer numtasks, rank, source, dest, tag, i,  ierr
    real*4 a(0:SIZE-1,0:SIZE-1), b(0:SIZE-1)
    integer stat(MPI_STATUS_SIZE)
    integer rowtype   ! required variable
    tag = 1

    ! Fortran stores this array in column major order
    data a  /1.0, 2.0, 3.0, 4.0, &
             5.0, 6.0, 7.0, 8.0,  &
             9.0, 10.0, 11.0, 12.0, &
             13.0, 14.0, 15.0, 16.0 /

    call MPI_INIT(ierr)
    call MPI_COMM_RANK(MPI_COMM_WORLD, rank, ierr)
    call MPI_COMM_SIZE(MPI_COMM_WORLD, numtasks, ierr)

    ! =====================================================================
    ! Step 1. Create an MPI Vector Type
    !    Summary:
    !      Make a new vector derived datatype.
    !
    !    Function Call:
    !      MPI_TYPE_VECTOR(INTEGER count,
    !                      INTEGER blocklength,
    !                      INTEGER stride,
    !                      INTEGER oldtype,
    !                      INTEGER newtype,
    !                      INTEGER ierr)
    !      
    !   Input Parameters:
    !       count
    !           number of blocks (nonnegative integer)
    !       blocklength
    !           number of elements in each block (nonnegative integer)
    !       stride
    !           number of elements between start of each block (integer)
    !       oldtype
    !           old datatype (handle)
    !       
    !    Output Parameters:
    !       newtype
    !          new datatype (handle)
    !
    ! TODO: create the vector data type
    
    ! TODO: commit the new derived datatype
    
    ! =====================================================================

    if (numtasks .eq. SIZE) then
        ! task 0 sends one element of rowtype to all tasks
        if (rank .eq. 0) then
            ! ===================================================================
            ! Step 2. Send vector data type using MPI_Send.
            !    Summary:
            !      Call MPI_Send 
            !
            !    Function Call:
            !      MPI_SEND(type(*) buf,
            !               INTEGER count,
            !               INTEGER datatype,
            !               INTEGER dest,
            !               INTEGER tag,
            !               INTEGER comm,
            !               INTEGER ierr)
            !      
            !   Input Parameters:
            !     buf
            !         initial address of send buffer (choice)
            !     count
            !         number of elements in send buffer (non-negative integer)
            !     datatype
            !         datatype of each send buffer element (handle)
            !     dest
            !         rank of destination (integer)
            !     tag
            !         message tag (integer)
            !     comm
            !         communicator (handle)
            !
            do i=0, numtasks-1
                ! TODO: send each COLUMN i of the array 'a' using the derived data type.
                
            end do
            ! ===================================================================
        endif

        ! all tasks receive rowtype data from task 0
        source = 0
        call MPI_RECV(b, SIZE, MPI_REAL, source, tag, MPI_COMM_WORLD, stat, ierr)
        print *, 'rank= ',rank,' b= ',b
    else
        print *, 'Must specify',SIZE,' processors.  Terminating.' 
    endif

    ! free datatype when done using it
    call MPI_TYPE_FREE(rowtype, ierr)
    call MPI_FINALIZE(ierr)

end program

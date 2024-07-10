program contiguous
    include 'mpif.h'

    integer SIZE
    parameter(SIZE=4)
    integer numtasks, rank, source, dest, tag, i,  ierr
    real*4 a(0:SIZE-1,0:SIZE-1), b(0:SIZE-1)
    integer stat(MPI_STATUS_SIZE)
    integer columntype   ! required variable
    tag = 1

    ! Fortran stores this array in column major order
    data a  /1.0, 2.0, 3.0, 4.0, &
             5.0, 6.0, 7.0, 8.0, &
             9.0, 10.0, 11.0, 12.0, & 
             13.0, 14.0, 15.0, 16.0 /

    call MPI_INIT(ierr)
    call MPI_COMM_RANK(MPI_COMM_WORLD, rank, ierr)
    call MPI_COMM_SIZE(MPI_COMM_WORLD, numtasks, ierr)

    ! =====================================================================
    ! Step 1. Create an MPI Contiguous Type
    !    Summary:
    !      Make a new contiguous derived datatype.
    !
    !    Function Call:
    !      MPI_TYPE_CONTIGUOUS(INTEGER count,
    !                          INTEGER oldtype,
    !                          INTEGER newtype,
    !                          INTEGER ierr)
    !      
    !   Input Parameters:
    !       count
    !           replication count (non-negative integer)
    !       oldtype
    !           old datatype (handle)
    !
    !   Output Parameters
    !       newtype
    !           new datatype (handle)
    !
    ! TODO: create the contiguous data type
    
    ! TODO: commit the new derived datatype 
    
    ! =====================================================================

    if (numtasks .eq. SIZE) then
        ! task 0 sends one element of columntype to all tasks
        if (rank .eq. 0) then
            ! ===================================================================
            ! Step 2. Send continguous data type using MPI_Send.
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
            do i=1, numtasks-1
                !  TODO: send each ROW i of the array 'a' using the derived data type.
            
            end do
            ! fill b for rank 0
            b = a(:,0)
            ! ===================================================================
        endif

        ! all tasks receive columntype data from task 0
        source = 0
        if(rank /= 0) call MPI_RECV(b, SIZE, MPI_REAL, source, tag, MPI_COMM_WORLD, stat, ierr)
        print *, 'rank= ',rank,' b= ',b
    else
        print *, 'Must specify',SIZE,' processors.  Terminating.' 
    endif

    ! free datatype when done using it
    call MPI_TYPE_FREE(columntype, ierr)
    call MPI_FINALIZE(ierr)

end program

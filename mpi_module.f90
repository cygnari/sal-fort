MODULE MPI_Module
    use mpi
    implicit none

    CONTAINS

    ! Array A is fully located at processor source, distribute so that each processor has a copy of A
    SUBROUTINE Distribute_Array(A, size, source, comm, p)
        use mpi
        integer, intent(in) :: size, source, comm, p
        integer :: win, ierr
        real(8), intent(inout) :: A(size)
        integer (kind=MPI_ADDRESS_KIND) :: lowerbound, sizeofreal, disp
        disp = 0

        CALL MPI_TYPE_GET_EXTENT(MPI_DOUBLE_PRECISION, lowerbound, sizeofreal, ierr)
        CALL MPI_WIN_Create(A, size*sizeofreal, MPI_DOUBLE_PRECISION, MPI_INFO_NULL, comm, win, ierr)
        CALL MPI_WIN_Fence(0, win, ierr)
        IF (p /= source) THEN
            CALL MPI_Get(A, size, MPI_DOUBLE_PRECISION, 0, disp, size, MPI_DOUBLE_PRECISION, win, ierr)
        END IF
        CALL MPI_WIN_Fence(0, win, ierr)
        CALL MPI_WIN_Free(win, ierr)
    END SUBROUTINE

    ! sum up array A that has entries across all processors
    SUBROUTINE Sum_Array(A, array_size, comm, id)
        use mpi
        integer, intent(in) :: array_size, comm, id
        integer :: win, ierr
        real(8), intent(inout) :: A(array_size)
        integer (kind=MPI_ADDRESS_KIND) :: lowerbound, sizeofreal, disp
        disp = 0

        CALL MPI_TYPE_GET_EXTENT(MPI_DOUBLE_PRECISION, lowerbound, sizeofreal, ierr)
        CALL MPI_WIN_Create(A, array_size*sizeofreal, MPI_DOUBLE_PRECISION, MPI_INFO_NULL, comm, win, ierr)
        CALL MPI_WIN_Fence(0, win, ierr)
        IF (id /= 0) THEN
            CALL MPI_Accumulate(A, array_size, MPI_DOUBLE_PRECISION, 0, disp, array_size, MPI_DOUBLE_PRECISION, MPI_SUM, win, ierr)
        END IF
        CALL MPI_WIN_Fence(0, win, ierr)
        IF (id /= 0) THEN
            CALL MPI_Get(A, array_size, MPI_DOUBLE_PRECISION, 0, disp, array_size, MPI_DOUBLE_PRECISION, win, ierr)
        END IF
        CALL MPI_WIN_Fence(0, win, ierr)
        CALL MPI_WIN_Free(win, ierr)
    END SUBROUTINE

    SUBROUTINE Gather_points(val, array, comm, p, id)
        ! each rank has val, construct array with vals from each processor
        use mpi
        integer, intent(in) :: val, comm, p, id
        integer, intent(out) :: array(:)
        integer :: ierr
        call MPI_Barrier(comm, ierr)
        call MPI_Allgather(val, 1, MPI_Integer, array, 1, MPI_Integer, comm, ierr)
        call MPI_Barrier(comm, ierr)
    END SUBROUTINE

    SUBROUTINE Gather_point_data(data, all_points, own_points, point_count, proc_points, starting_points, comm, p, id)
        use mpi
        real(8), intent(in) :: data(:)
        real(8), intent(out) :: all_points(:)
        integer, intent(in) :: starting_points(:), proc_points(:)
        integer, intent(in) :: comm, p, id, own_points, point_count
        integer :: ierr
        call MPI_Barrier(comm, ierr)
        call MPI_Allgatherv(data, own_points, MPI_DOUBLE_PRECISION, all_points, proc_points, starting_points, &
                            MPI_DOUBLE_PRECISION, comm, ierr)
        call MPI_Barrier(comm, ierr)
    END SUBROUTINE

    ! do precomputation for ssh precomputation
    SUBROUTINE Communicate_sshs_pre(unowned_points, own_points, starting_points_procs, &
                                comm, P, ID, points_to_give_proc, needed_points_proc, points_send_proc)
        use mpi
        integer, intent(in) :: unowned_points(:,:), own_points, comm, P, ID
        integer, intent(in) :: starting_points_procs(:), needed_points_proc(:)
        integer :: ierr, count, i, j, request, maxpoints, shift
        integer, allocatable, intent(out) :: points_to_give_proc(:)
        integer, allocatable, intent(out) :: points_send_proc(:,:)
        integer, allocatable :: send_reqs(:), recv_reqs(:), status(:,:)
        call MPI_Barrier(comm, ierr)
        ! first have each processor communicate how many points it needs from each other processor
        allocate(points_to_give_proc(P), source=-1)
        allocate(send_reqs(P))
        allocate(recv_reqs(P))
        allocate(status(MPI_STATUS_SIZE,P))
        
        ! print *, 'here 1'
        DO i = 1, P ! each processor gathers the number of points each processor needs 
            call MPI_Gather(needed_points_proc(i), 1, MPI_INTEGER, points_to_give_proc, 1, MPI_INTEGER, &
                                i-1, comm, ierr)
            ! rank i-1 needs needed_points_proc(i) points from rank ID
            ! rank ID gives points_to_give_proc(i) to rank i-1
        END DO
        call MPI_Barrier(comm, ierr)
        ! print *, 'here 2'
        maxpoints = 0
        DO i = 1, P
            maxpoints = max(maxpoints, points_to_give_proc(i))
        END DO
        allocate(points_send_proc(maxpoints, P), source=-1)
        ! each processor gathers the points that all the other processors need
        DO i = 1, P ! each processor sends the points it needs to rank i-1
            call MPI_Isend(unowned_points(:,i), needed_points_proc(i), MPI_Integer, i-1, i-1, comm, send_reqs(i), ierr)
        END DO
        ! print *, 'here 3'
        DO i = 1, P ! each processor receives the point indices it was sent from rank i-1
            call MPI_Irecv(points_send_proc(:,i), points_to_give_proc(i), MPI_Integer, i-1, ID, comm, recv_reqs(i), ierr)
        END DO
        ! print *, 'here 4'
        call MPI_Waitall(P, send_reqs, status, ierr)
        call MPI_Waitall(P, recv_reqs, status, ierr)
        ! print *, 'here 5'
        call MPI_Barrier(comm, ierr)
        ! shift points to local indices
        shift = starting_points_procs(ID+1)-1
        DO i = 1, P
            DO j = 1, points_to_give_proc(i)
                points_send_proc(j, i) = points_send_proc(j, i) - shift
            END DO
        END DO
    END SUBROUTINE

    ! communicate sshs
    SUBROUTINE Communicate_sshs(sshs, e_sshs, points_send_proc, points_to_give_proc, needed_points_proc, &
                                unowned_points, comm, P, ID)
        real(8), intent(in) :: sshs(:)
        real(8), intent(out) :: e_sshs(:)
        integer, intent(in) :: points_send_proc(:,:), points_to_give_proc(:), needed_points_proc(:), comm, P, ID
        integer, intent(in) :: unowned_points
        real(8), allocatable :: sshs_to_send(:,:), sshs_to_receive(:,:)
        integer :: maxpoints, i, j, ierr, index
        integer, allocatable :: send_reqs(:), recv_reqs(:), status(:,:)
        allocate(send_reqs(P))
        allocate(recv_reqs(P))
        allocate(status(MPI_STATUS_SIZE,P))
        maxpoints = 0
        DO i = 1, P
            maxpoints = max(maxpoints, points_to_give_proc(i))
        END DO
        allocate(sshs_to_send(maxpoints, P), source=0.0_8)
        maxpoints = 0
        DO i = 1, P
            maxpoints = max(maxpoints, needed_points_proc(i))
        END DO
        allocate(sshs_to_receive(maxpoints, P), source=0.0_8)
        DO i = 1, P
            DO j = 1, points_to_give_proc(i)
                ! copy the sshs that proc i needs
                sshs_to_send(j,i) = sshs(points_send_proc(j,i))
            END DO
        END DO
        ! send sshs
        DO i = 1, P
            call MPI_Isend(sshs_to_send(:,i), points_to_give_proc(i), MPI_DOUBLE_PRECISION, i-1, i-1, comm, send_reqs(i), ierr)
        END DO
        ! recv sshs
        DO i = 1, P
            call MPI_Irecv(sshs_to_receive(:,i), needed_points_proc(i), MPI_DOUBLE_PRECISION, i-1, ID, comm, recv_reqs(i), ierr)
        END DO
        call MPI_Waitall(P, send_reqs, status, ierr)
        call MPI_Waitall(P, recv_reqs, status, ierr)
        call MPI_Barrier(comm, ierr)
        index = 0
        DO i = 1, P
            DO j = 1, needed_points_proc(i)
                index = index + 1
                e_sshs(index) = sshs_to_receive(j, i)
            END DO
        END DO
        ! print *, ID, index
    END SUBROUTINE
END MODULE
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
    SUBROUTINE Sum_Array(A, array_size, comm, p)
        use mpi
        integer, intent(in) :: array_size, comm, p
        integer :: win, ierr
        real(8), intent(inout) :: A(array_size)
        integer (kind=MPI_ADDRESS_KIND) :: lowerbound, sizeofreal, disp
        disp = 0

        CALL MPI_TYPE_GET_EXTENT(MPI_DOUBLE_PRECISION, lowerbound, sizeofreal, ierr)
        CALL MPI_WIN_Create(A, array_size*sizeofreal, MPI_DOUBLE_PRECISION, MPI_INFO_NULL, comm, win, ierr)
        CALL MPI_WIN_Fence(0, win, ierr)
        IF (p /= 0) THEN
            CALL MPI_Accumulate(A, array_size, MPI_DOUBLE_PRECISION, 0, disp, array_size, MPI_DOUBLE_PRECISION, MPI_SUM, win, ierr)
        END IF
        CALL MPI_WIN_Fence(0, win, ierr)
        IF (p /= 0) THEN
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
END MODULE
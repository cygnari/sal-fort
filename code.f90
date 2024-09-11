program main
    use mpi
    use tree_module
    use cube_sphere_module
    implicit none

    real(8) :: theta, area, t1, t2
    integer :: point_count, ierr, process_rank, total_ranks, i, cluster_thresh, interp_degree
    real(8), allocatable :: areas(:), xs(:), ys(:), zs(:), sshs(:), sal(:), sal_x(:), sal_y(:)
    type(cube_panel), allocatable :: tree_panels(:)
    type(interaction_pair), allocatable :: interaction_list(:)
    NAMELIST /params/ theta, point_count, cluster_thresh, interp_degree

    call MPI_INIT(ierr)
    call MPI_COMM_SIZE(MPI_COMM_WORLD, total_ranks, ierr)
    call MPI_COMM_RANK(MPI_COMM_WORLD, process_rank, ierr)

    IF(cluster_thresh == -1) THEN
        cluster_thresh = interp_degree * interp_degree + 2
    END IF

    IF(process_rank == 0) THEN
        CALL CPU_TIME(t1)
    END IF

    open (file='namelist.nml', unit=9)
    read (nml=params, unit=9)
    close (9)

    allocate(areas(point_count)) ! initialize arrays of area, x co, y co, z co, ssh
    allocate(xs(point_count))
    allocate(ys(point_count))
    allocate(zs(point_count))
    allocate(sshs(point_count))

    IF (process_rank == 0) THEN ! only read in values on one processor
        open (file='./data-files/2313486_mpas_grid_areas.csv', unit=10) ! read in values from data files
        open (file='./data-files/2313486_mpas_grid_x.csv', unit=11)
        open (file='./data-files/2313486_mpas_grid_y.csv', unit=12)
        open (file='./data-files/2313486_mpas_grid_z.csv', unit=13)
        open (file='./data-files/2313486_mpas_grid_ssh.csv', unit=14)

        DO i = 1, point_count
            read (10, *) areas(i)
            read (11, *) xs(i)
            read (12, *) ys(i)
            read (13, *) zs(i)
            read (14, *) sshs(i)
        END DO

        close (10)
        close (11)
        close (12)
        close (13)
        close (14)
    END IF

    ! distribute data values to every other processor
    call Distribute_Array(areas, point_count, 0, MPI_COMM_WORLD, process_rank)
    call Distribute_Array(xs, point_count, 0, MPI_COMM_WORLD, process_rank)
    call Distribute_Array(ys, point_count, 0, MPI_COMM_WORLD, process_rank)
    call Distribute_Array(zs, point_count, 0, MPI_COMM_WORLD, process_rank)
    call Distribute_Array(sshs, point_count, 0, MPI_COMM_WORLD, process_rank)

    IF (process_rank == 0) THEN
        CALL CPU_TIME(t2)
        print *, 'initialize points time : ',(t2 - t1)
        CALL CPU_TIME(t1)
    END IF

    ! construct the cubed sphere tree
    call tree_traverse(tree_panels, xs, ys, zs, point_count, cluster_thresh, process_rank)

    IF (process_rank == 0) THEN
        CALL CPU_TIME(t2)
        print *, 'tree traverse time : ',(t2 - t1)
        CALL CPU_TIME(t1)
    END IF

    print *, 'total cube panels: ', size(tree_panels)

    ! perform dual tree traversal to compute interaction lists

    call dual_tree_traversal(interaction_list, tree_panels, theta, cluster_thresh)

    IF (process_rank == 0) THEN
        CALL CPU_TIME(t2)
        print *, 'DTT time: ',(t2 - t1)
        CALL CPU_TIME(t1)
    END IF

    print *, 'total interactions: ', size(interaction_list)

    ! fast summation to approximate the convolution

    IF (process_rank == 0) THEN
        CALL CPU_TIME(t2)
        print *, 'Integration time: ',(t2 - t1)
        CALL CPU_TIME(t1)
    END IF

    call MPI_FINALIZE(ierr)

end program

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

SUBROUTINE Sum_Array(A, size, comm, p)
    use mpi
    integer, intent(in) :: size, comm, p
    integer :: win, ierr
    real(8), intent(inout) :: A(size)
    integer (kind=MPI_ADDRESS_KIND) :: lowerbound, sizeofreal, disp
    disp = 0

    CALL MPI_TYPE_GET_EXTENT(MPI_DOUBLE_PRECISION, lowerbound, sizeofreal, ierr)
    CALL MPI_WIN_Create(A, size*sizeofreal, MPI_DOUBLE_PRECISION, MPI_INFO_NULL, comm, win, ierr)
    CALL MPI_WIN_Fence(0, win, ierr)
    CALL MPI_Accumulate(A, size, MPI_DOUBLE_PRECISION, 0, disp, size, MPI_DOUBLE_PRECISION, MPI_SUM, win, ierr)
    CALL MPI_WIN_Fence(0, win, ierr)
    IF (p /= 0) THEN
        CALL MPI_Get(A, size, MPI_DOUBLE_PRECISION, 0, disp, size, MPI_DOUBLE_PRECISION, win, ierr)
    END IF
    CALL MPI_WIN_Fence(0, win, ierr)
    CALL MPI_WIN_Free(win, ierr)
END SUBROUTINE

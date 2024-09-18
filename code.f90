program main
    use mpi
    use tree_module
    use cube_sphere_module
    use fast_sum_module
    use mpi_module
    implicit none

    real(8) :: theta, area, t1, t2
    integer :: point_count, ierr, process_rank, total_ranks, i, cluster_thresh, interp_degree
    real(8), allocatable :: areas(:), xs(:), ys(:), zs(:), sshs(:), sal(:), sal_x(:), sal_y(:)
    type(cube_panel), allocatable :: tree_panels(:)
    type(interaction_pair), allocatable :: interaction_list(:), own_interactions(:)
    NAMELIST /params/ theta, point_count, cluster_thresh, interp_degree
    ! theta is fast sum theta parameter, cluster_thresh is point count threshold for cluster
    ! interp_degree is the degree of the barycentric lagrange interpolation to do

    call MPI_INIT(ierr)
    call MPI_COMM_SIZE(MPI_COMM_WORLD, total_ranks, ierr)
    call MPI_COMM_RANK(MPI_COMM_WORLD, process_rank, ierr)

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
        open(file='./data-files/2313486_mpas_grid_areas.csv', unit=10) ! read in values from data files
        open(file='./data-files/2313486_mpas_grid_x.csv', unit=11)
        open(file='./data-files/2313486_mpas_grid_y.csv', unit=12)
        open(file='./data-files/2313486_mpas_grid_z.csv', unit=13)
        open(file='./data-files/2313486_mpas_grid_ssh.csv', unit=14)

        DO i = 1, point_count
            read(10, *) areas(i)
            read(11, *) xs(i)
            read(12, *) ys(i)
            read(13, *) zs(i)
            read(14, *) sshs(i)
        END DO

        close(10)
        close(11)
        close(12)
        close(13)
        close(14)
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

    ! construct the cubed sphere tree structure in tree_panels
    call tree_traverse(tree_panels, xs, ys, zs, point_count, cluster_thresh, process_rank)

    IF (process_rank == 0) THEN
        CALL CPU_TIME(t2)
        print *, 'tree traverse time : ',(t2 - t1)
        CALL CPU_TIME(t1)
        print *, 'total cube panels: ', size(tree_panels)
    END IF

    ! perform dual tree traversal to compute interaction lists
    call dual_tree_traversal(interaction_list, tree_panels, theta, cluster_thresh)

    IF (process_rank == 0) THEN
        CALL CPU_TIME(t2)
        print *, 'DTT time: ',(t2 - t1)
        CALL CPU_TIME(t1)
        print *, 'total interactions: ', size(interaction_list)
    END IF

    ! fast summation to approximate the convolution
    allocate(sal_x(point_count), source=0.0_8)
    allocate(sal_y(point_count), source=0.0_8)

    CALL fast_sum(sal_x, sal_y, interaction_list, tree_panels, xs, ys, zs, areas, sshs, point_count, & 
                    interp_degree, process_rank, total_ranks)

    ! replace these two with reproducing sums if needed
    CALL Sum_Array(sal_x, point_count, MPI_COMM_WORLD, process_rank)
    CALL Sum_array(sal_y, point_count, MPI_COMM_WORLD, process_rank)

    IF (process_rank == 0) THEN
        CALL CPU_TIME(t2)
        print *, 'Integration time: ',(t2 - t1)
        CALL CPU_TIME(t1)

        open(file='./run-output/2313486_sal_x.csv', status='replace', unit = 15)
        open(file='./run-output/2313486_sal_y.csv', status='replace', unit = 16)
        DO i = 1, point_count
            write(15, *) sal_x(i)
            write(16, *) sal_y(i)
        END DO
        close(15)
        close(16)

        CALL CPU_TIME(t2)
        print *, 'Output write time: ', (t2-t1)
    END IF

    call MPI_FINALIZE(ierr)

end program
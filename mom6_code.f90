program main
    use mpi
    use tree_module
    use cube_sphere_module
    use fast_sum_module
    use mpi_module
    implicit none

    real(8) :: theta, area, t1, t2, min_x, max_x, min_y, max_y, min_z, max_z, x, y, z, a, ssh
    integer :: point_count, ierr, process_rank, total_ranks, i, cluster_thresh, interp_degree, own_points
    integer :: max_level, level, j, unowned_needed_points, duplicates
    real(8), allocatable :: areas(:), xs(:), ys(:), zs(:), sshs(:), sal(:), sal_x(:), sal_y(:), e_areas(:)
    real(8), allocatable :: areas_t(:), xs_t(:), ys_t(:), zs_t(:), sshs_t(:), e_xs(:), e_ys(:), e_zs(:), e_sshs(:)
    integer, allocatable :: points_per_proc(:), starting_point_proc(:), point_proc_id(:), points_panels(:,:)
    integer, allocatable :: needed_points_proc(:), unowned_points(:,:)
    type(cube_panel), allocatable :: tree_panels(:)
    type(interaction_pair), allocatable :: interaction_list(:), own_interactions(:)
    NAMELIST /params/ theta, point_count, cluster_thresh, interp_degree
    ! theta is fast sum theta parameter, cluster_thresh is point count threshold for cluster
    ! interp_degree is the degree of the barycentric lagrange interpolation to do

    call MPI_INIT(ierr)
    call MPI_COMM_SIZE(MPI_COMM_WORLD, total_ranks, ierr)
    call MPI_COMM_RANK(MPI_COMM_WORLD, process_rank, ierr)

    ! domain decomposition
    IF (total_ranks == 1) THEN
        min_x = -2.0_8
        max_x = 2.0_8
        min_y = -2.0_8
        max_y = 2.0_8
        min_z = -2.0_8
        max_z = 2.0_8
    ELSE IF (total_ranks == 2) THEN
        ! split x coordinate 
        min_x = -1.0_8+ibits(process_rank, 0, 1)*1.0_8
        max_x = 0.0_8+ibits(process_rank, 0, 1)*1.0_8
        min_y = -2.0_8
        max_y = 2.0_8
        min_z = -2.0_8
        max_z = 2.0_8
    ELSE IF (total_ranks == 4) THEN
        ! split x and y coordinates
        min_x = -1.0_8+ibits(process_rank, 0, 1)*1.0_8
        max_x = 0.0_8+ibits(process_rank, 0, 1)*1.0_8
        min_y = -1.0_8+ibits(process_rank, 1, 1)*1.0_8
        max_y = 0.0_8+ibits(process_rank, 1, 1)*1.0_8
        min_z = -2.0_8
        max_z = 2.0_8
    ELSE IF (total_ranks == 8) THEN
        ! split x y z coordinates
        min_x = -1.0_8+ibits(process_rank, 0, 1)*1.0_8
        max_x = 0.0_8+ibits(process_rank, 0, 1)*1.0_8
        min_y = -1.0_8+ibits(process_rank, 1, 1)*1.0_8
        max_y = 0.0_8+ibits(process_rank, 1, 1)*1.0_8
        min_z = -1.0_8+ibits(process_rank, 2, 1)*1.0_8
        max_z = 0.0_8+ibits(process_rank, 2, 1)*1.0_8
    END IF

    IF(process_rank == 0) THEN
        CALL CPU_TIME(t1)
    END IF

    open (file='namelist.nml', unit=9)
    read (nml=params, unit=9)
    close (9)

    IF (process_rank == 0) THEN
        print *, 'theta: ', theta, ' point count: ', point_count, ' cluster_thresh: ', cluster_thresh, &
                ' interp_degree: ', interp_degree
    END IF

    allocate(areas_t(point_count), source=0.0_8) ! initialize arrays of area, x co, y co, z co, ssh
    allocate(xs_t(point_count), source=0.0_8)
    allocate(ys_t(point_count), source=0.0_8)
    allocate(zs_t(point_count), source=0.0_8)
    allocate(sshs_t(point_count), source=0.0_8)

    own_points = 0

    ! each processor reads its own points, processors don't start with all the points
    open(file='./data-files/2313486_mpas_grid_areas.csv', unit=10) ! read in values from data files
    open(file='./data-files/2313486_mpas_grid_x.csv', unit=11)
    open(file='./data-files/2313486_mpas_grid_y.csv', unit=12)
    open(file='./data-files/2313486_mpas_grid_z.csv', unit=13)
    open(file='./data-files/2313486_mpas_grid_ssh.csv', unit=14)

    DO i = 1, point_count
        read(10, *) a
        read(11, *) x
        read(12, *) y
        read(13, *) z
        read(14, *) ssh
        ! print *, x, y, z, a, ssh
        IF ((x >= min_x) .and. (x < max_x) .and. (y >= min_y) .and. (y < max_y) .and. (z >= min_z) .and. (z < max_z)) THEN
            ! print *, 'here'
            own_points = own_points + 1
            xs_t(own_points) = x
            ys_t(own_points) = y
            zs_t(own_points) = z
            areas_t(own_points) = a
            sshs_t(own_points) = ssh
        END IF
    END DO

    close(10)
    close(11)
    close(12)
    close(13)
    close(14)

    allocate(areas(own_points))
    allocate(xs(own_points))
    allocate(ys(own_points))
    allocate(zs(own_points))
    allocate(sshs(own_points))

    DO i = 1, own_points
        areas(i) = areas_t(i)
        areas_t(i) = 0.0_8
        xs(i) = xs_t(i)
        xs_t(i) = 0.0_8
        ys(i) = ys_t(i)
        ys_t(i) = 0.0_8
        zs(i) = zs_t(i)
        zs_t(i) = 0.0_8
        sshs(i) = sshs_t(i)
        sshs_t(i) = 0.0_8
    END DO

    IF (process_rank == 0) THEN
        CALL CPU_TIME(t2)
        print *, 'read points time : ',(t2 - t1)
        CALL CPU_TIME(t1)
    END IF

    call MPI_Barrier(MPI_COMM_WORLD, ierr)

    print *, 'rank: ', process_rank, ' points: ', own_points

    ! construct tree of source points

    allocate(points_per_proc(total_ranks))
    allocate(starting_point_proc(total_ranks))
    allocate(point_proc_id(point_count))

    call Gather_points(own_points, points_per_proc, MPI_COMM_WORLD, total_ranks, process_rank)

    starting_point_proc(1) = 1
    DO i = 2, total_ranks
        starting_point_proc(i) = starting_point_proc(i-1)+points_per_proc(i-1)
    END DO

    call Gather_point_data(xs, xs_t, own_points, point_count, points_per_proc, starting_point_proc, &
                            MPI_COMM_WORLD, total_ranks, process_rank)
    call Gather_point_data(ys, ys_t, own_points, point_count, points_per_proc, starting_point_proc, &
                            MPI_COMM_WORLD, total_ranks, process_rank)
    call Gather_point_data(zs, zs_t, own_points, point_count, points_per_proc, starting_point_proc, &
                            MPI_COMM_WORLD, total_ranks, process_rank)

    ! construct cubed sphere tree
    call tree_traverse(tree_panels, xs_t, ys_t, zs_t, point_count, cluster_thresh, process_rank)

    IF (process_rank == 0) THEN
        CALL CPU_TIME(t2)
        print *, 'tree construction time : ',(t2 - t1)
        CALL CPU_TIME(t1)
        print *, 'total cube panels: ', size(tree_panels)
    END IF

    ! determine which panels own points are in

    max_level = tree_panels(size(tree_panels))%level
    allocate(points_panels(max_level, own_points), source = -1)

    call assign_points_to_panels(xs, ys, zs, points_panels, tree_panels, own_points)

    IF (process_rank == 0) THEN
        CALL CPU_TIME(t2)
        print *, 'point panel assignment time : ',(t2 - t1)
        CALL CPU_TIME(t1)
    END IF

    ! perform single tree traversal to compute interaction lists

    call single_tree_traversal(interaction_list, tree_panels, theta, cluster_thresh, xs, ys, zs, own_points)

    IF (process_rank == 0) THEN
        CALL CPU_TIME(t2)
        print *, 'interaction list computation time: ',(t2 - t1)
        CALL CPU_TIME(t1)
    END IF
    call MPI_Barrier(MPI_COMM_WORLD, ierr)
    print *, 'total interactions: ', size(interaction_list), ' average interactions per point: ', &
                1.0_8*size(interaction_list)/own_points

    ! compute source points from other processors that are needed

    call unowned_sources(starting_point_proc, interaction_list, tree_panels, total_ranks, process_rank, unowned_points, &
                            needed_points_proc, point_count)

    unowned_needed_points = 0
    DO i = 1, total_ranks
        ! print *, process_rank, i, needed_points_proc(i)
        unowned_needed_points = unowned_needed_points + needed_points_proc(i)
    END DO
    print *, 'process ', process_rank, 'needs ', unowned_needed_points, ' unowned points'
    call MPI_Barrier(MPI_COMM_WORLD, ierr)
    ! IF (process_rank == 3) THEN
    !     DO i = 1, total_ranks
    !         DO j = 1, needed_points_proc(i)
    !             print *, i, j, unowned_points(j, i)
    !         END DO
    !     END DO
    ! END IF

    allocate(e_xs(unowned_needed_points))
    allocate(e_ys(unowned_needed_points))
    allocate(e_zs(unowned_needed_points))
    allocate(e_sshs(unowned_needed_points))
    allocate(e_areas(unowned_needed_points))

    

    ! ! fast summation to approximate the convolution
    ! allocate(sal_x(point_count), source=0.0_8)
    ! allocate(sal_y(point_count), source=0.0_8)

    ! CALL fast_sum(sal_x, sal_y, interaction_list, tree_panels, xs, ys, zs, areas, sshs, point_count, & 
    !                 interp_degree, process_rank, total_ranks)

    ! ! replace these two with reproducing sums if needed
    ! CALL Sum_Array(sal_x, point_count, MPI_COMM_WORLD, process_rank)
    ! CALL Sum_array(sal_y, point_count, MPI_COMM_WORLD, process_rank)

    ! IF (process_rank == 0) THEN
    !     CALL CPU_TIME(t2)
    !     print *, 'Integration time: ',(t2 - t1)
    !     CALL CPU_TIME(t1)

    !     open(file='./run-output/2313486_sal_x.csv', status='replace', unit = 15)
    !     open(file='./run-output/2313486_sal_y.csv', status='replace', unit = 16)
    !     DO i = 1, point_count
    !         write(15, *) sal_x(i)
    !         write(16, *) sal_y(i)
    !     END DO
    !     close(15)
    !     close(16)

    !     CALL CPU_TIME(t2)
    !     print *, 'Output write time: ', (t2-t1)
    ! END IF

    call MPI_FINALIZE(ierr)

end program
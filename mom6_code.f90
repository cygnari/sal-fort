program main
    use mpi
    use tree_module
    use cube_sphere_module
    use mom6_fast_sum_module
    use mpi_module
    implicit none

    real(8) :: theta, area, t1, t2, min_x, max_x, min_y, max_y, min_z, max_z, x, y, z, a, ssh
    integer :: point_count, ierr, process_rank, total_ranks, i, cluster_thresh, interp_degree, own_points
    integer :: max_level, level, j, unowned_needed_points, duplicates, interact_count, k
    real(8), allocatable :: areas(:), xs(:), ys(:), zs(:), sshs(:), sal(:), sal_x(:), sal_y(:), e_areas(:)
    ! xs, ys, zs are the locations of the points that a processor owns
    real(8), allocatable :: areas_t(:), xs_t(:), ys_t(:), zs_t(:), sshs_t(:), e_xs(:), e_ys(:), e_zs(:), e_sshs(:)
    ! xs_t, ys_t, zs_t are all the points across all the ranks, used only for initial tree traversal/interaction list computation
    ! e_xs, e_ys, e_zs are the points that a processor doesn't own but needs for pp interactions
    ! e_sshs are the corresponding unowned areas that need to be communicated
    real(8), allocatable :: proxy_source_weights(:), sal_xs(:), sal_ys(:)
    integer, allocatable :: points_per_proc(:), starting_point_proc(:), point_proc_id(:), points_panels(:,:)
    ! points_panels(:,i) is the panels containing point i
    integer, allocatable :: needed_points_proc(:), unowned_points(:,:), points_to_give_proc(:)
    ! unowned points are the points needed for PP interactions that a processor does not own
    ! 2d array, second dim is processor index
    ! needed_points_proc is how many points are needed from each other rank
    ! points_to_give_proc is how many points are to be given to each other rank
    integer, allocatable :: points_needed_from_proc(:), points_send_proc(:,:)
    ! points_send_proc is a list of sshs to send each processor
    ! 2d array, second dim is processor index
    type(cube_panel), allocatable :: tree_panels(:)
    type(interaction_pair), allocatable :: pc_interaction_list(:), pp_interaction_list(:)
    NAMELIST /params/ theta, point_count, cluster_thresh, interp_degree
    ! theta is fast sum theta parameter, cluster_thresh is point count threshold for cluster
    ! interp_degree is the degree of the barycentric lagrange interpolation to do

    call MPI_INIT(ierr)
    call MPI_COMM_SIZE(MPI_COMM_WORLD, total_ranks, ierr)
    call MPI_COMM_RANK(MPI_COMM_WORLD, process_rank, ierr)

    ! domain decomposition, relatively naive way of splitting up the earth
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
    call Gather_point_data(areas, areas_t, own_points, point_count, points_per_proc, starting_point_proc, &
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
    allocate(points_panels(max_level+1, own_points), source = -1)

    call assign_points_to_panels(xs, ys, zs, points_panels, tree_panels, own_points)

    IF (process_rank == 0) THEN
        CALL CPU_TIME(t2)
        print *, 'point panel assignment time : ',(t2 - t1)
        CALL CPU_TIME(t1)
    END IF

    ! perform single tree traversal to compute interaction lists

    call single_tree_traversal(pc_interaction_list, pp_interaction_list, tree_panels, theta, cluster_thresh, &
                                xs, ys, zs, own_points)

    IF (process_rank == 0) THEN
        CALL CPU_TIME(t2)
        print *, 'interaction list computation time: ',(t2 - t1)
        CALL CPU_TIME(t1)
    END IF
    call MPI_Barrier(MPI_COMM_WORLD, ierr)

    interact_count = size(pc_interaction_list) + size(pp_interaction_list)

    print *, 'total interactions: ', interact_count, ' average interactions per point: ', &
                1.0_8*interact_count/own_points
    print *, 'pc interactions: ', size(pc_interaction_list), ' pp interactions: ', size(pp_interaction_list)

    ! compute source points from other processors that are needed

    call unowned_sources(starting_point_proc, pp_interaction_list, tree_panels, total_ranks, process_rank, unowned_points, &
                            needed_points_proc, point_count)

    unowned_needed_points = 0
    DO i = 1, total_ranks
        unowned_needed_points = unowned_needed_points + needed_points_proc(i)
    END DO
    print *, 'process ', process_rank, 'needs ', unowned_needed_points, ' unowned points'
    call MPI_Barrier(MPI_COMM_WORLD, ierr)

    allocate(e_xs(unowned_needed_points), source=0.0_8)
    allocate(e_ys(unowned_needed_points), source=0.0_8)
    allocate(e_zs(unowned_needed_points), source=0.0_8)
    allocate(e_sshs(unowned_needed_points), source=0.0_8)
    allocate(e_areas(unowned_needed_points), source=0.0_8)

    call relabel_copy_unowned(e_xs, e_ys, e_zs, e_areas, unowned_points, xs_t, ys_t, zs_t, areas_t, tree_panels, &
                                    total_ranks, process_rank, unowned_needed_points, own_points, &
                                    starting_point_proc)

    IF (process_rank .ne. 0) THEN
        deallocate(xs_t)
        deallocate(ys_t)
        deallocate(zs_t)
        deallocate(areas_t)
        deallocate(sshs_t)
    END IF

    IF (process_rank == 0) THEN
        CALL CPU_TIME(t2)
        print *, 'relabel time: ',(t2 - t1)
        CALL CPU_TIME(t1)
    END IF

    ! communicate needed SSHs

    call Communicate_sshs_pre(unowned_points, own_points, starting_point_proc, MPI_COMM_WORLD, &
                                total_ranks, process_rank, points_to_give_proc, needed_points_proc, &
                                points_send_proc)

    call MPI_Barrier(MPI_COMM_WORLD, ierr)

    call Communicate_sshs(sshs, e_sshs, points_send_proc, points_to_give_proc, needed_points_proc, &
                                unowned_needed_points, MPI_COMM_WORLD, total_ranks, process_rank)

    call MPI_Barrier(MPI_COMM_WORLD, ierr)
    IF (process_rank == 0) THEN
        CALL CPU_TIME(t2)
        print *, 'sshs communication time: ',(t2 - t1)
        CALL CPU_TIME(t1)
    END IF

    ! fast summation to approximate the convolution
    ! first compute proxy source weights
    allocate(proxy_source_weights((interp_degree+1)*(interp_degree+1)*size(tree_panels)), source=0.0_8)

    call proxy_source_compute(proxy_source_weights, tree_panels, points_panels, xs, ys, zs, sshs, areas, interp_degree, &
                                    own_points)

    ! IF (process_rank == 0) THEN
    !     DO i = 1, size(proxy_source_weights)
    !         print *, proxy_source_weights(i)
    !     END DO
    ! END IF

    call Sum_Array(proxy_source_weights, (interp_degree+1)*(interp_degree+1)*size(tree_panels), MPI_COMM_WORLD, process_rank)

    call MPI_Barrier(MPI_COMM_WORLD, ierr)

    ! IF (process_rank == 0) THEN
    !     DO i = 1, 4
    !         print *, proxy_source_weights(i)
    !     END DO
    ! END IF

    IF (process_rank == 0) THEN
        CALL CPU_TIME(t2)
        print *, 'proxy source computation time: ',(t2 - t1)
        CALL CPU_TIME(t1)
    END IF

    allocate(sal_x(own_points), source=0.0_8)
    allocate(sal_y(own_points), source=0.0_8)

    call pc_compute(sal_x, sal_y, pc_interaction_list, tree_panels, xs, ys, zs, interp_degree, proxy_source_weights)

    call MPI_Barrier(MPI_COMM_WORLD, ierr)
    IF (process_rank == 0) THEN
        CALL CPU_TIME(t2)
        print *, 'pc computation time: ',(t2 - t1)
        CALL CPU_TIME(t1)
    END IF

    call pp_compute(sal_x, sal_y, pp_interaction_list, tree_panels, xs, ys, zs, areas, sshs, e_xs, e_ys, e_zs, e_areas, &
                    e_sshs, own_points)

    call MPI_Barrier(MPI_COMM_WORLD, ierr)
    IF (process_rank == 0) THEN
        CALL CPU_TIME(t2)
        print *, 'pp computation time: ',(t2 - t1)
        CALL CPU_TIME(t1)
    END IF

    allocate(sal_xs(point_count), source=0.0_8)
    allocate(sal_ys(point_count), source=0.0_8)

    ! print *, 'here 1', sal_y(1)

    call Gather_point_data(sal_x, sal_xs, own_points, point_count, points_per_proc, starting_point_proc, &
                            MPI_COMM_WORLD, total_ranks, process_rank)
    call Gather_point_data(sal_y, sal_ys, own_points, point_count, points_per_proc, starting_point_proc, &
                            MPI_COMM_WORLD, total_ranks, process_rank)

    ! print *, 'here 2', sal_ys(1)

    call MPI_Barrier(MPI_COMM_WORLD, ierr)
    IF (process_rank == 0) THEN
        CALL CPU_TIME(t2)
        print *, 'Communication time: ',(t2 - t1)
        CALL CPU_TIME(t1)

        open(file='./run-output/2313486_sal_x.csv', status='replace', unit = 15)
        open(file='./run-output/2313486_sal_y.csv', status='replace', unit = 16)
        open(file='./run-output/2313486_x.csv', status='replace', unit = 17)
        open(file='./run-output/2313486_y.csv', status='replace', unit = 18)
        open(file='./run-output/2313486_z.csv', status='replace', unit = 19)
        open(file='./run-output/2313486_areas.csv', status='replace', unit = 20)
        DO i = 1, point_count
            write(15, *) sal_xs(i)
            write(16, *) sal_ys(i)
            write(17, *) xs_t(i)
            write(18, *) ys_t(i)
            write(19, *) zs_t(i)
            write(20, *) areas_t(i)
        END DO
        close(15)
        close(16)

        CALL CPU_TIME(t2)
        print *, 'Output write time: ', (t2-t1)
    END IF

    call MPI_FINALIZE(ierr)

end program
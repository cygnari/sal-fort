MODULE Tree_Module
    use cube_sphere_module
    implicit none

    type :: cube_panel ! cubed sphere tree panel face data type
        integer :: level = 0 ! level 0 is the base level
        logical :: is_leaf = .true. ! all leaves start out as leaves, become not leaf when refined
        integer :: id
        integer :: parent_panel = -1
        integer :: child_panel_1 = -1
        integer :: child_panel_2 = -1
        integer :: child_panel_3 = -1
        integer :: child_panel_4 = -1
        integer :: face 
        real(8) :: min_xi ! xi and eta are the angle coordinates on the face of the cube
        real(8) :: mid_xi ! based on the equiangular gnomonic cubed sphere
        real(8) :: max_xi
        real(8) :: min_eta
        real(8) :: mid_eta
        real(8) :: max_eta
        real(8) :: radius ! distance from center of panel to corner
        integer, allocatable :: points_inside(:)
        integer, allocatable :: relabeled_points_inside(:)
        integer :: panel_point_count = 0

        contains
            procedure :: contains_point
    end type cube_panel

    type :: interaction_pair ! data type to represent a P/C-P/C interaction between two cube panels
        integer :: index_target
        integer :: index_source 
        integer :: interact_type ! 0 for PP, 1 for PC, 2 for CP, 3 for CC
    end type interaction_pair

    CONTAINS 

    ! type bound procedure, determine if panel self contains point x, y, z
    logical function contains_point(self, x, y, z) result(contains)
        class(cube_panel), intent(in) :: self
        real(8), intent(in) :: x, y, z
        integer :: face
        real(8) :: xi, eta
        face = face_from_xyz(x, y, z)
        IF (face == self%face) THEN
            call xieta_from_xyz(x, y, z, xi, eta, face)
            IF ((xi >= self%min_xi) .and. (xi < self%max_xi) .and. (eta >= self%min_eta) &
                            .and. (eta < self%max_eta)) THEN
                contains = .true.
            ELSE
                contains = .false.
            END IF
        ELSE 
            contains = .false.
        END IF
    end function

    ! performs tree traversal of the points (xs, ys, zs), giving an array tree_panels of cube_panels
    SUBROUTINE tree_traverse(tree_panels, xs, ys, zs, point_count, cluster_thresh, process_rank)
        type(cube_panel), allocatable, intent(out) :: tree_panels(:)
        type(cube_panel), allocatable :: tree_panels_temp(:)
        integer, intent(in) :: point_count, cluster_thresh, process_rank
        real(8), intent(in) :: xs(point_count), ys(point_count), zs(point_count)
        integer :: i, face, total, j, panel_count, index, which_panel, count, k
        real(8) :: pi, xval, yval, zval, cube_xi, cube_eta, min_xi, max_xi, mid_xi, min_eta, mid_eta, max_eta, xi, eta
        real(8) :: x1, x2, x3, y1, y2, y3, d1, d2, d3, d4
        logical :: found, point_found
        integer, allocatable :: curr_loc(:), temp(:)

        pi = 4.D0*DATAN(1.D0)

        allocate(tree_panels_temp(max(point_count,6)))
        allocate(curr_loc(6))

        ! initialize the six top level cube panels
        DO i = 1, 6
            tree_panels_temp(i)%id = i
            tree_panels_temp(i)%face = i
            tree_panels_temp(i)%min_xi = -pi/4.D0
            tree_panels_temp(i)%mid_xi = 0.D0
            tree_panels_temp(i)%max_xi = pi/4.D0
            tree_panels_temp(i)%min_eta = -pi/4.D0
            tree_panels_temp(i)%mid_eta = 0.D0
            tree_panels_temp(i)%max_eta = pi/4.D0
            allocate(tree_panels_temp(i)%points_inside(point_count))
            tree_panels_temp(i)%panel_point_count = 0
            curr_loc(i) = 0
        END DO

        panel_count = 6

        ! assign points to the top level cube panels
        DO i = 1, point_count
            xval = xs(i)
            yval = ys(i)
            zval = zs(i)
            face = face_from_xyz(xval, yval, zval)
            curr_loc(face) = curr_loc(face) + 1
            tree_panels_temp(face)%panel_point_count = tree_panels_temp(face)%panel_point_count + 1
            tree_panels_temp(face)%points_inside(curr_loc(face)) = i
        END DO

        ! sanity check to make sure that the points were assigned correctly to the six top level panels
        total = 0
        DO i = 1, 6
            total = total + tree_panels_temp(i)%panel_point_count
        END DO

        IF (process_rank == 0) THEN
            IF (total /= point_count) THEN
                print *, 'Total assigned points ', total, ' does not equal point count ', point_count
            END IF
        END IF

        ! resize the points_inside arrays of the top level panels
        DO i = 1, 6
            allocate(temp(curr_loc(i)))
            DO j = 1, tree_panels_temp(i)%panel_point_count
                temp(j) = tree_panels_temp(i)%points_inside(j)
            END DO
            call move_alloc(from=temp, to=tree_panels_temp(i)%points_inside)
            curr_loc(i) = 0
        END DO

        ! iterate through the panels, if a panel is a cluster, refine it and create four new clusters
        i = 1
        DO WHILE (i <= panel_count)
            ! first check if panel needs to be divided
            count = tree_panels_temp(i)%panel_point_count
            IF ((count >= cluster_thresh) .and. (tree_panels_temp(i)%is_leaf)) THEN
                ! panel is a leaf and has many points => refine
                tree_panels_temp(i)%is_leaf = .false.
                min_xi = tree_panels_temp(i)%min_xi
                mid_xi = tree_panels_temp(i)%mid_xi
                max_xi = tree_panels_temp(i)%max_xi
                min_eta = tree_panels_temp(i)%min_eta
                mid_eta = tree_panels_temp(i)%mid_eta
                max_eta = tree_panels_temp(i)%max_eta
                ! four new panels are index panel_count+1 through panel_count+4
                DO j = 1, 4
                    index = panel_count+j
                    ! print *, index, i
                    tree_panels_temp(panel_count+j)%parent_panel = i
                    tree_panels_temp(panel_count+j)%level = tree_panels_temp(i)%level + 1
                    tree_panels_temp(panel_count+j)%face = tree_panels_temp(i)%face
                    tree_panels_temp(panel_count+j)%id = panel_count+j
                    allocate(tree_panels_temp(index)%points_inside(count))
                END DO
                ! coordinates of four new panels in angle space
                tree_panels_temp(panel_count+1)%min_xi = min_xi; tree_panels_temp(panel_count+1)%min_eta = min_eta
                tree_panels_temp(panel_count+1)%max_xi = mid_xi; tree_panels_temp(panel_count+1)%max_eta = mid_eta
                tree_panels_temp(panel_count+1)%mid_xi = 0.5*(min_xi+mid_xi)
                tree_panels_temp(panel_count+1)%mid_eta = 0.5*(min_eta+mid_eta)

                tree_panels_temp(panel_count+2)%min_xi = mid_xi; tree_panels_temp(panel_count+2)%min_eta = min_eta
                tree_panels_temp(panel_count+2)%max_xi = max_xi; tree_panels_temp(panel_count+2)%max_eta = mid_eta
                tree_panels_temp(panel_count+2)%mid_xi = 0.5*(mid_xi+max_xi)
                tree_panels_temp(panel_count+2)%mid_eta = 0.5*(min_eta+mid_eta)

                tree_panels_temp(panel_count+3)%min_xi = mid_xi; tree_panels_temp(panel_count+3)%min_eta = mid_eta
                tree_panels_temp(panel_count+3)%max_xi = max_xi; tree_panels_temp(panel_count+3)%max_eta = max_eta
                tree_panels_temp(panel_count+3)%mid_xi = 0.5*(mid_xi+max_xi)
                tree_panels_temp(panel_count+3)%mid_eta = 0.5*(mid_eta+max_eta)

                tree_panels_temp(panel_count+4)%min_xi = min_xi; tree_panels_temp(panel_count+4)%min_eta = mid_eta
                tree_panels_temp(panel_count+4)%max_xi = mid_xi; tree_panels_temp(panel_count+4)%max_eta = max_eta
                tree_panels_temp(panel_count+4)%mid_xi = 0.5*(min_xi+mid_xi)
                tree_panels_temp(panel_count+4)%mid_eta = 0.5*(mid_eta+max_eta)

                tree_panels_temp(i)%child_panel_1 = panel_count+1
                tree_panels_temp(i)%child_panel_2 = panel_count+2
                tree_panels_temp(i)%child_panel_3 = panel_count+3
                tree_panels_temp(i)%child_panel_4 = panel_count+4

                DO j = 1, tree_panels_temp(i)%panel_point_count
                    ! loop through points contained in parent panel, assign to sub panels
                    index = tree_panels_temp(i)%points_inside(j)
                    xval = xs(index); yval = ys(index); zval = zs(index)
                    call xieta_from_xyz(xval, yval, zval, xi, eta, tree_panels_temp(i)%face)
                    IF (xi < mid_xi) THEN
                        IF (eta < mid_eta) THEN
                            which_panel = 1
                        ELSE 
                            which_panel = 4
                        END IF
                    ELSE 
                        IF (eta < mid_eta) THEN
                            which_panel = 2
                        ELSE 
                            which_panel = 3
                        END IF
                    END IF
                    curr_loc(which_panel) = curr_loc(which_panel) + 1
                    tree_panels_temp(panel_count+which_panel)%panel_point_count = & 
                                                tree_panels_temp(panel_count+which_panel)%panel_point_count + 1
                    tree_panels_temp(panel_count+which_panel)%points_inside(curr_loc(which_panel)) = index
                END DO

                ! resize points_inside array
                DO j = 1, 4
                    allocate(temp(curr_loc(j)))
                    DO k = 1, tree_panels_temp(panel_count+j)%panel_point_count
                        temp(k) = tree_panels_temp(panel_count+j)%points_inside(k)
                    END DO
                    call move_alloc(from=temp, to=tree_panels_temp(panel_count+j)%points_inside)
                    curr_loc(j) = 0
                END DO

                ! sanity check to make sure all points are assigned
                total = tree_panels_temp(panel_count+1)%panel_point_count
                total = total + tree_panels_temp(panel_count+2)%panel_point_count
                total = total + tree_panels_temp(panel_count+3)%panel_point_count
                total = total + tree_panels_temp(panel_count+4)%panel_point_count
                IF (total /= tree_panels_temp(i)%panel_point_count) THEN
                    print *, 'Error in refining tree at parent panel ', i
                END IF
                panel_count = panel_count + 4
            END IF
            i = i + 1
        END DO

        ! tree_panels_temp is the wrong size so move everything over to the correctly sized tree_panels
        allocate(tree_panels(panel_count))
        DO i = 1, panel_count
            tree_panels(i) = tree_panels_temp(i)
            tree_panels(i)%relabeled_points_inside = tree_panels(i)%points_inside
            min_xi = tree_panels(i)%min_xi
            mid_xi = tree_panels(i)%mid_xi
            max_xi = tree_panels(i)%max_xi
            min_eta = tree_panels(i)%min_eta
            mid_eta = tree_panels(i)%mid_eta
            max_eta = tree_panels(i)%max_eta
            ! compute the furthest distance from panel center to vertex for each panel
            call xyz_from_xieta(x1, x2, x3, mid_xi, mid_eta, tree_panels(i)%face)
            call xyz_from_xieta(y1, y2, y3, min_xi, min_eta, tree_panels(i)%face)
            d1 = ACOS(MAX(MIN(x1*y1+x2*y2+x3*y3, 1.0_8), -1.0_8))
            call xyz_from_xieta(y1, y2, y3, min_xi, max_eta, tree_panels(i)%face)
            d2 = ACOS(MAX(MIN(x1*y1+x2*y2+x3*y3, 1.0_8), -1.0_8))
            call xyz_from_xieta(y1, y2, y3, max_xi, max_eta, tree_panels(i)%face)
            d3 = ACOS(MAX(MIN(x1*y1+x2*y2+x3*y3, 1.0_8), -1.0_8))
            call xyz_from_xieta(y1, y2, y3, max_xi, min_eta, tree_panels(i)%face)
            d4 = ACOS(MAX(MIN(x1*y1+x2*y2+x3*y3, 1.0_8), -1.0_8))
            tree_panels(i)%radius = MAX(d1, d2, d3, d4)
        END DO

    END SUBROUTINE

    ! find which panels each point is in
    SUBROUTINE assign_points_to_panels(xs, ys, zs, points_panels, tree_panels, own_points)
        real(8), intent(in) :: xs(:), ys(:), zs(:)
        integer, intent(out) :: points_panels(:,:)
        type(cube_panel), intent(in) :: tree_panels(:)
        integer, intent(in) :: own_points
        integer :: level, i, j, face
        real(8) :: xi, eta, xi1, xi2, eta1, eta2

        DO i = 1, own_points
            level = 1
            j = 1
            jloop: DO ! do loop over tree panels
                IF (j == -1) THEN
                    exit jloop
                ELSE IF (tree_panels(j)%contains_point(xs(i), ys(i), zs(i))) THEN
                    ! point i is contained in panel j
                    points_panels(level, i) = j
                    level = level + 1
                    j = tree_panels(j)%child_panel_1
                ELSE
                    j = j + 1
                END IF
            END DO jloop
        END DO
    END SUBROUTINE

    ! dual tree traversal to compute interaction list
    ! consider the tree of targets and tree of sources, work one panel at a time and determine if interact or refine
    SUBROUTINE DUAL_TREE_TRAVERSAL(interaction_list, tree_panels, param_theta, param_cluster_thresh)
        type(interaction_pair), allocatable, intent(out) :: interaction_list(:)
        type(cube_panel), intent(in) :: tree_panels(:)
        real(8), intent(in) :: param_theta
        integer, intent(in) :: param_cluster_thresh
        integer :: interaction_count, curr_loc, i, j, i_t, i_s, c_t, c_s, tree_traverse_count, int_type
        real(8) :: x1t, x2t, x3t, x1s, x2s, x3s, dist, separation
        integer, allocatable :: target_index(:), source_index(:)
        type(interaction_pair), allocatable :: interaction_lists_temp(:)

        ! initialize arrays, experience tells us that there should be around 32*panel_count interactions
        allocate(target_index(size(tree_panels)*64))
        allocate(source_index(size(tree_panels)*64))
        allocate(interaction_lists_temp(size(tree_panels)*64))
        tree_traverse_count = 0

        DO i = 1, 6
            DO j = 1, 6
                tree_traverse_count = tree_traverse_count + 1
                target_index(tree_traverse_count) = i
                source_index(tree_traverse_count) = j
            END DO
        END DO
        
        curr_loc = 1
        interaction_count = 0
        DO WHILE (curr_loc <= tree_traverse_count)
            i_t = target_index(curr_loc)
            i_s = source_index(curr_loc)
            c_t = tree_panels(i_t)%panel_point_count
            c_s = tree_panels(i_s)%panel_point_count
            IF ((c_t > 0) .and. (c_s > 0)) THEN ! check both panels have points
                call xyz_from_xieta(x1t, x2t, x3t, tree_panels(i_t)%mid_xi, tree_panels(i_t)%mid_eta, tree_panels(i_t)%face)
                call xyz_from_xieta(x1s, x2s, x3s, tree_panels(i_s)%mid_xi, tree_panels(i_s)%mid_eta, tree_panels(i_s)%face)
                dist = ACOS(MIN(MAX(x1t*x1s+x2t*x2s+x3t*x3s, -1.0_8), 1.0_8))
                separation = (tree_panels(i_t)%radius+tree_panels(i_s)%radius)/dist
                IF ((dist > 0) .and. (separation < param_theta)) THEN
                    ! two panels are well separated
                    int_type = 0
                    interaction_count = interaction_count + 1
                    IF (c_t > param_cluster_thresh) THEN
                        int_type = int_type + 2
                    END IF
                    IF (c_s > param_cluster_thresh) THEN
                        int_type = int_type + 1
                    END IF
                    interaction_lists_temp(interaction_count)%index_target = i_t
                    interaction_lists_temp(interaction_count)%index_source = i_s
                    interaction_lists_temp(interaction_count)%interact_type = int_type
                ELSE
                    ! two panels are not well separated
                    IF ((c_t < param_cluster_thresh) .and. (c_s < param_cluster_thresh)) THEN
                        ! both panels have few points
                        interaction_count = interaction_count + 1
                        interaction_lists_temp(interaction_count)%index_target = i_t
                        interaction_lists_temp(interaction_count)%index_source = i_s
                        interaction_lists_temp(interaction_count)%interact_type = 0
                    ELSE IF (tree_panels(i_t)%is_leaf .and. tree_panels(i_s)%is_leaf) THEN
                        ! both panels are leaves, cannot be refined
                        interaction_count = interaction_count + 1
                        interaction_lists_temp(interaction_count)%index_target = i_t
                        interaction_lists_temp(interaction_count)%index_source = i_s
                        interaction_lists_temp(interaction_count)%interact_type = 0
                    ELSE IF (tree_panels(i_t)%is_leaf) THEN
                        ! target panel is leaf, refine source
                        target_index(tree_traverse_count+1) = i_t
                        target_index(tree_traverse_count+2) = i_t
                        target_index(tree_traverse_count+3) = i_t
                        target_index(tree_traverse_count+4) = i_t
                        source_index(tree_traverse_count+1) = tree_panels(i_s)%child_panel_1
                        source_index(tree_traverse_count+2) = tree_panels(i_s)%child_panel_2
                        source_index(tree_traverse_count+3) = tree_panels(i_s)%child_panel_3
                        source_index(tree_traverse_count+4) = tree_panels(i_s)%child_panel_4
                        tree_traverse_count = tree_traverse_count + 4
                    ELSE IF (tree_panels(i_s)%is_leaf) THEN
                        ! source panel is leaf, refine target
                        source_index(tree_traverse_count+1) = i_s
                        source_index(tree_traverse_count+2) = i_s
                        source_index(tree_traverse_count+3) = i_s
                        source_index(tree_traverse_count+4) = i_s
                        target_index(tree_traverse_count+1) = tree_panels(i_t)%child_panel_1
                        target_index(tree_traverse_count+2) = tree_panels(i_t)%child_panel_2
                        target_index(tree_traverse_count+3) = tree_panels(i_t)%child_panel_3
                        target_index(tree_traverse_count+4) = tree_panels(i_t)%child_panel_4
                        tree_traverse_count = tree_traverse_count + 4
                    ELSE 
                        ! neither source nor target is a leaf, refine the panel with more points
                        IF (c_t >= c_s) THEN
                            ! target panel has more points, refine target
                            source_index(tree_traverse_count+1) = i_s
                            source_index(tree_traverse_count+2) = i_s
                            source_index(tree_traverse_count+3) = i_s
                            source_index(tree_traverse_count+4) = i_s
                            target_index(tree_traverse_count+1) = tree_panels(i_t)%child_panel_1
                            target_index(tree_traverse_count+2) = tree_panels(i_t)%child_panel_2
                            target_index(tree_traverse_count+3) = tree_panels(i_t)%child_panel_3
                            target_index(tree_traverse_count+4) = tree_panels(i_t)%child_panel_4
                            tree_traverse_count = tree_traverse_count + 4
                        ELSE 
                            ! source panel has more points, refine source
                            target_index(tree_traverse_count+1) = i_t
                            target_index(tree_traverse_count+2) = i_t
                            target_index(tree_traverse_count+3) = i_t
                            target_index(tree_traverse_count+4) = i_t
                            source_index(tree_traverse_count+1) = tree_panels(i_s)%child_panel_1
                            source_index(tree_traverse_count+2) = tree_panels(i_s)%child_panel_2
                            source_index(tree_traverse_count+3) = tree_panels(i_s)%child_panel_3
                            source_index(tree_traverse_count+4) = tree_panels(i_s)%child_panel_4
                            tree_traverse_count = tree_traverse_count + 4
                        END IF
                    END IF
                END IF
            END IF
            curr_loc = curr_loc + 1
        END DO

        ! copy interactions from temp list 
        allocate(interaction_list(interaction_count))
        DO i = 1, interaction_count
            interaction_list(i) = interaction_lists_temp(i)
        END DO
    END SUBROUTINE

    ! performs tree traversal only with the source tree for individual target particles
    SUBROUTINE SINGLE_TREE_TRAVERSAL(pc_interaction_list, pp_interaction_list, tree_panels, param_theta, param_cluster_thresh, &
                                        xs, ys, zs, own_points)
        type(interaction_pair), allocatable, intent(out) :: pc_interaction_list(:), pp_interaction_list(:)
        type(cube_panel), intent(in) :: tree_panels(:)
        real(8), intent(in) :: param_theta, xs(:), ys(:), zs(:)
        integer, intent(in) :: param_cluster_thresh, own_points
        integer :: interaction_count, curr_loc, i, j, i_t, i_s, c_t, c_s, tree_traverse_count, int_type, pc_count, pp_count, k, l
        real(8) :: x1t, x2t, x3t, x1s, x2s, x3s, dist, separation
        integer, allocatable :: source_index(:)
        type(interaction_pair), allocatable :: interaction_lists_temp(:)
        logical :: well_separated

        allocate(source_index(size(tree_panels)))
        allocate(interaction_lists_temp(own_points*128))

        interaction_count = 0
        pp_count = 0
        pc_count = 0
        DO i = 1, own_points ! loop over own target points
            tree_traverse_count = 6
            DO j = 1, 6
                source_index(j) = j
            END DO
            curr_loc = 1
            x1t = xs(i)
            x2t = ys(i)
            x3t = zs(i)
            DO WHILE (curr_loc <= tree_traverse_count) ! go through the source panels
                i_s = source_index(curr_loc)
                c_s = tree_panels(i_s)%panel_point_count
                IF (c_s > 0) THEN ! check that source panel has points
                    ! two ways of checking separation
                    ! way 1, compute theta from distance between center of panel and target
                    call xyz_from_xieta(x1s, x2s, x3s, tree_panels(i_s)%mid_xi, tree_panels(i_s)%mid_eta, tree_panels(i_s)%face)
                    dist = ACOS(MIN(MAX(x1t*x1s+x2t*x2s+x3t*x3s, -1.0_8), 1.0_8))
                    separation = tree_panels(i_s)%radius/dist
                    well_separated = .false.
                    IF ((dist > 0) .and. (separation < param_theta)) THEN
                        well_separated = .true.
                    END IF
                    ! way 2, check if target i is in panel i_s
                    ! well_separated = .not. tree_panels(i_s)%contains_point(x1t, x2t, x3t)
                    IF (well_separated) THEN
                        ! well separated, do cluster interaction
                        interaction_count = interaction_count + 1
                        interaction_lists_temp(interaction_count)%index_target = i
                        interaction_lists_temp(interaction_count)%index_source = i_s
                        interaction_lists_temp(interaction_count)%interact_type = 1
                        pc_count = pc_count + 1
                        ! print *, 1, i, j
                    ELSE IF (tree_panels(i_s)%is_leaf) THEN
                        ! source panel is a leaf, cannot refine further
                        interaction_count = interaction_count + 1
                        interaction_lists_temp(interaction_count)%index_target = i
                        interaction_lists_temp(interaction_count)%index_source = i_s
                        interaction_lists_temp(interaction_count)%interact_type = 0
                        pp_count = pp_count + 1
                        ! print *, 0, i, j
                    ELSE 
                        ! source panel is not a leaf and is not well separated, refine 
                        source_index(tree_traverse_count+1) = tree_panels(i_s)%child_panel_1
                        source_index(tree_traverse_count+2) = tree_panels(i_s)%child_panel_2
                        source_index(tree_traverse_count+3) = tree_panels(i_s)%child_panel_3
                        source_index(tree_traverse_count+4) = tree_panels(i_s)%child_panel_4
                        tree_traverse_count = tree_traverse_count + 4
                    END IF
                END IF
                curr_loc = curr_loc + 1
            END DO
        END DO

        ! print *, pc_count, pp_count

        allocate(pc_interaction_list(pc_count))
        allocate(pp_interaction_list(pp_count))
        k = 0
        l = 0
        DO i = 1, interaction_count
            IF (interaction_lists_temp(i)%interact_type == 1) THEN
                k = k + 1
                pc_interaction_list(k) = interaction_lists_temp(i)
            ELSE
                l = l + 1
                pp_interaction_list(l) = interaction_lists_temp(i)
            END IF
        END DO
    END SUBROUTINE

    ! see which points are needed for pp interactions that a processor does not own
    SUBROUTINE unowned_sources(starting_point_proc, pp_interaction_list, tree_panels, P, ID, unowned_points, &
                                points_per_proc, total_points)
        integer, intent(in) :: starting_point_proc(:), P, ID, total_points
        type(interaction_pair), intent(in) :: pp_interaction_list(:)
        type(cube_panel), intent(in) :: tree_panels(:)
        integer, intent(out), allocatable :: unowned_points(:,:), points_per_proc(:)
        integer, allocatable :: unowned_points_temp(:,:)
        integer :: points_needed, i, own_min, own_max, j, i_s, i_sp, k, kmax, kmin, l
        logical :: found

        points_needed = 0
        allocate(points_per_proc(P), source=0)
        allocate(unowned_points_temp(total_points, P))

        own_min = starting_point_proc(ID+1)
        IF (ID+1==P) THEN
            ! last proc
            own_max = total_points
        ELSE
            own_max = starting_point_proc(ID+2)-1
        END IF

        iloop: DO i = 1, size(pp_interaction_list)
            ! particle particle interaction
            i_s = pp_interaction_list(i)%index_source
            jloop: DO j = 1, tree_panels(i_s)%panel_point_count
                ! loop over points in source panel
                i_sp = tree_panels(i_s)%points_inside(j)
                IF ((i_sp < own_min) .or. (i_sp > own_max)) THEN
                    ! point is not owned
                    points_needed = points_needed + 1
                    kloop: DO k = 1, P 
                        ! loop over points, find which processor owns it
                        kmin = starting_point_proc(k)
                        IF (k == P) THEN
                            kmax = total_points
                        ELSE
                            kmax = starting_point_proc(k+1)-1
                        END IF
                        IF ((i_sp <= kmax) .and. (i_sp >= kmin)) THEN
                            ! point i_sp owned by processor k (MPI rank k-1)
                            found = .false.
                            lloop: DO l = 1, points_per_proc(k)
                                ! check if l is already contained in unowned_points_temp
                                IF (unowned_points_temp(l, k) == i_sp) THEN
                                    found = .true.
                                END IF
                            END DO lloop
                            IF (.not. found) THEN
                                points_per_proc(k) = points_per_proc(k) + 1
                                unowned_points_temp(points_per_proc(k), k) = i_sp
                            END IF
                            exit kloop
                        END IF
                    END DO kloop
                END IF
            END DO jloop
        END DO iloop

        kmax = 0
        DO i = 1, P
            kmax = max(kmax, points_per_proc(i))
        END DO

        allocate(unowned_points(kmax, P), source=-1)
        DO j = 1, P
            DO i = 1, points_per_proc(j)
                unowned_points(i, j) = unowned_points_temp(i, j)
            END DO
        END DO

    END SUBROUTINE

    ! for all the unowned source points that a processor needs for particle particle interactions
    ! relabel them so that they are contiguous in e_xs/e_ys/e_zs using values from xs/ys/zs_t 
    ! renumber them in tree_panels as well
    SUBROUTINE relabel_copy_unowned(e_xs, e_ys, e_zs, e_areas, unowned_points, xs_t, ys_t, zs_t, areas_t, tree_panels, &
                                    P, ID, unowned_points_count, own_points, starting_points)
        real(8), intent(out) :: e_xs(unowned_points_count), e_ys(unowned_points_count), e_zs(unowned_points_count), &
                                e_areas(unowned_points_count)
        integer, intent(in) :: unowned_points(:,:), P, ID, unowned_points_count, own_points, starting_points(:)
        real(8), intent(in) :: xs_t(:), ys_t(:), zs_t(:), areas_t(:)
        type(cube_panel), intent(inout) :: tree_panels(:)
        ! integer, intent(out) :: relabel_map(:)
        integer :: curr_index, i, j, i_p, k, l, min_i, max_i
        logical :: found

        min_i = starting_points(ID+1)
        max_i = min_i+own_points

        DO i = 1, size(tree_panels)
            DO j = 1, tree_panels(i)%panel_point_count
                i_p = tree_panels(i)%points_inside(j)
                IF ((i_p >= min_i) .and. (i_p < max_i)) THEN
                    tree_panels(i)%relabeled_points_inside(j) = i_p-min_i+1
                END IF
            END DO
        END DO

        curr_index = 0
        DO i = 1, P
            DO j = 1, size(unowned_points, 1)
                ! loop over points in unowned arrays
                ! print *, i, j
                i_p = unowned_points(j, i)
                IF (i_p .ne. -1) THEN
                    curr_index = curr_index + 1
                    e_xs(curr_index) = xs_t(i_p)
                    e_ys(curr_index) = ys_t(i_p)
                    e_zs(curr_index) = zs_t(i_p)
                    e_areas(curr_index) = areas_t(i_p)
                    k = 1
                    ! relabel_map(i_p) = curr_index
                    treeloop: DO ! do loop over tree panels
                        IF (k == -1) THEN ! leaf index
                            exit treeloop
                        ELSE
                            found = .false.
                            lloop: DO l = 1, tree_panels(k)%panel_point_count
                                IF (i_p == tree_panels(k)%points_inside(l)) THEN
                                    found = .true.
                                    tree_panels(k)%relabeled_points_inside(l) = own_points+curr_index
                                    exit lloop
                                END IF
                            END DO lloop
                            IF (found) THEN
                                k = tree_panels(k)%child_panel_1
                            ELSE
                                k = k + 1
                            END IF
                        END IF
                    END DO treeloop
                END IF  
            END DO
        END DO
    END SUBROUTINE

END MODULE Tree_Module
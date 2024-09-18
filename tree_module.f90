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
        integer :: panel_point_count = 0
    end type cube_panel

    type :: interaction_pair ! data type to represent a P/C-P/C interaction between two cube panels
        integer :: index_target
        integer :: index_source 
        integer :: interact_type ! 0 for PP, 1 for PC, 2 for CP, 3 for CC
    end type interaction_pair

    CONTAINS 

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

        allocate(tree_panels_temp(point_count))
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

    ! dual tree traversal to compute interaction list
    ! consider the tree of targets and tree of sources, work one panel at a time and determine if interact or refine
    SUBROUTINE DUAL_TREE_TRAVERSAL(interaction_list, tree_panels, param_theta, param_cluster_thresh)
        type(interaction_pair), allocatable, intent(out) :: interaction_list(:)
        type(cube_panel), intent(in) :: tree_panels(:)
        real(8), intent(in) :: param_theta
        integer, intent(in) :: param_cluster_thresh
        integer :: interaction_count, curr_loc, i, j, i_t, i_s, c_t, c_s, tree_traverse_count, int_type
        real(8) :: x1t, x2t, x3t, x1s, x2s, x3s, dist, separation
        integer, allocatable :: target_index(:), source_index(:), interaction_type(:)
        type(interaction_pair), allocatable :: interaction_lists_temp(:)

        ! initialize arrays, experience tells us that there should be around 32*panel_count interactions
        allocate(target_index(size(tree_panels)*64))
        allocate(source_index(size(tree_panels)*64))
        allocate(interaction_type(size(tree_panels)*64))
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

END MODULE Tree_Module
MODULE FAST_SUM_MODULE
    use tree_module
    use bli_module
    implicit none

    CONTAINS

    SUBROUTINE sal_grad_gfunc(tx, ty, tz, sx, sy, sz, sal_x, sal_y)
        real(8), intent(in) :: tx, ty, tz, sx, sy, sz
        real(8), intent(out) :: sal_x, sal_y
        real(8) :: g, mp, sqrtp, pi, cons, sqp, p1, p2, x32, val, mp2

        pi = 4.D0*DATAN(1.D0)
        ! cons = -3.0*1035.0/(3*6371000*4.0*pi*5517.0)
        cons = -2.343256857908601e-9_8
        ! cons = -7.029770573725803e-9_8

        sal_x = 0.0_8
        sal_y = 0.0_8
        IF ((ABS(tz - 1.0_8) > 1e-15_8) .and. (ABS(tz+1.0_8) > 1e-15_8)) THEN
            g = MAX(MIN(tx*sx+ty*sy+tz*sz, 1.0_8), -1.0_8) ! floating point check
            mp = 2.0_8-2.0_8*g
            sqp = SQRT(mp)
            p1 = (1.0_8-6.21196_8)/(sqp*mp+1e-16_8)
            p2 = (2.7_8+6_8)*(2_8*g+sqp) / (2_8*(g*g-1.0_8)+1e-16_8)
            val = (p1+p2)*cons
            x32 = tz*tz
            mp2 = SQRT(1.0_8-x32)
            sal_x = (sz*(1-x32)-tz*(tx*sx+ty*sy))/mp2*val
            sal_y = (tx*sy-ty*sx)/mp2*val
        END IF
    END SUBROUTINE

    SUBROUTINE pp_interact(sal_longrad, sal_latgrad, target_panel, source_panel, xs, ys, zs, area, ssh, interp_degree)
        real(8), intent(in) :: xs(:), ys(:), zs(:), area(:), ssh(:)
        type(cube_panel), intent(in) :: target_panel, source_panel
        integer, intent(in) :: interp_degree
        real(8), intent(inout) :: sal_longrad(:), sal_latgrad(:)
        real(8), allocatable :: sxs(:), sys(:), szs(:), sshs(:), areas(:)
        integer :: target_i, source_j, source_count, i, j
        real(8) :: tx, ty, tz, sx, sy, sz, sal_x, sal_y

        source_count = source_panel%panel_point_count
        allocate(sxs(source_count))
        allocate(sys(source_count))
        allocate(szs(source_count))
        allocate(sshs(source_count))
        allocate(areas(source_count))

        DO j= 1, source_count
            source_j = source_panel%points_inside(j)
            sxs(j) = xs(source_j)
            sys(j) = ys(source_j)
            szs(j) = zs(source_j)
            sshs(j) = ssh(source_j)
            areas(j) = area(source_j)
        END DO

        DO i = 1, target_panel%panel_point_count
            target_i = target_panel%points_inside(i)
            tx = xs(target_i)
            ty = ys(target_i)
            tz = zs(target_i)
            DO j = 1, source_count
                sx = sxs(j)
                sy = sys(j)
                sz = szs(j)
                call sal_grad_gfunc(tx, ty, tz, sx, sy, sz, sal_x, sal_y)
                sal_longrad(target_i) = sal_longrad(target_i) + sal_x*areas(j)*sshs(j)
                sal_latgrad(target_i) = sal_latgrad(target_i) + sal_y*areas(j)*sshs(j)
            END DO
        END DO
    END SUBROUTINE

    SUBROUTINE pc_interact(sal_longrad, sal_latgrad, target_panel, source_panel, xs, ys, zs, area, ssh, interp_degree)
        real(8), intent(in) :: xs(:), ys(:), zs(:), area(:), ssh(:)
        type(cube_panel), intent(in) :: target_panel, source_panel
        integer, intent(in) :: interp_degree
        real(8), intent(inout) :: sal_longrad(:), sal_latgrad(:)
    END SUBROUTINE

    SUBROUTINE cp_interact(sal_longrad, sal_latgrad, target_panel, source_panel, xs, ys, zs, area, ssh, interp_degree)
        real(8), intent(in) :: xs(:), ys(:), zs(:), area(:), ssh(:)
        type(cube_panel), intent(in) :: target_panel, source_panel
        integer, intent(in) :: interp_degree
        real(8), intent(inout) :: sal_longrad(:), sal_latgrad(:)
    END SUBROUTINE

    SUBROUTINE cc_interact(sal_longrad, sal_latgrad, target_panel, source_panel, xs, ys, zs, area, ssh, interp_degree)
        real(8), intent(in) :: xs(:), ys(:), zs(:), area(:), ssh(:)
        type(cube_panel), intent(in) :: target_panel, source_panel
        integer, intent(in) :: interp_degree
        real(8), intent(inout) :: sal_longrad(:), sal_latgrad(:)
    END SUBROUTINE

    SUBROUTINE fast_sum(sal_longrad, sal_latgrad, interaction_list, cube_panels, xs, ys, zs, area, ssh, point_count, &
                        interp_degree, rank, process_count)
        real(8), intent(in) :: xs(:), ys(:), zs(:), area(:), ssh(:)
        type(cube_panel), intent(in) :: cube_panels(:)
        type(interaction_pair), intent(in) :: interaction_list(:)
        integer, intent(in) :: point_count, interp_degree, rank, process_count
        real(8), intent(inout) :: sal_longrad(:), sal_latgrad(:)
        integer :: interaction_count, i_t, i_s, i, type, interactions

        interaction_count = size(interaction_list)
        interactions = 0

        DO i = 1, interaction_count
            IF (MOD(i, process_count) == rank) THEN
                i_t = interaction_list(i)%index_target
                i_s = interaction_list(i)%index_source
                type = interaction_list(i)%interact_type
                interactions = interactions + 1
                IF (type == 0) THEN
                    call pp_interact(sal_longrad, sal_latgrad, cube_panels(i_t), cube_panels(i_s), xs, ys, zs, &
                                        area, ssh, interp_degree)
                ELSE IF (type == 1) THEN
                    call pc_interact(sal_longrad, sal_latgrad, cube_panels(i_t), cube_panels(i_s), xs, ys, zs, &
                                        area, ssh, interp_degree)
                ELSE IF (type == 2) THEN
                    call cp_interact(sal_longrad, sal_latgrad, cube_panels(i_t), cube_panels(i_s), xs, ys, zs, &
                                        area, ssh, interp_degree)
                ELSE IF (type == 3) THEN
                    call cc_interact(sal_longrad, sal_latgrad, cube_panels(i_t), cube_panels(i_s), xs, ys, zs, &
                                        area, ssh, interp_degree)
                END IF
            END IF
        END DO

        print *, 'Rank: ', rank, ' interactions: ', interactions
    END SUBROUTINE
END MODULE
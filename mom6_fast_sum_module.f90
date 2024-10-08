MODULE MOM6_FAST_SUM_MODULE
    use tree_module
    use bli_module
    use cube_sphere_module
    implicit none

    CONTAINS

    SUBROUTINE sal_grad_gfunc(tx, ty, tz, sx, sy, sz, sal_x, sal_y)
        real(8), intent(in) :: tx, ty, tz, sx, sy, sz
        real(8), intent(out) :: sal_x, sal_y
        real(8) :: g, mp, sqrtp, cons, sqp, p1, p2, x32, val, mp2

        cons = -2.343256857908601e-9_8
        ! cons = -7.029770573725803e-9_8

        sal_x = 0.0_8
        sal_y = 0.0_8
        IF ((ABS(tz - 1.0_8) > 1e-15_8) .and. (ABS(tz+1.0_8) > 1e-15_8)) THEN
            g = MAX(MIN(tx*sx+ty*sy+tz*sz, 1.0_8), -1.0_8) ! floating point check
            mp = 2.0_8-2.0_8*g
            sqp = SQRT(mp)
            p1 = (1.0_8-6.21196_8)/(sqp*mp+1e-16_8)
            p2 = (2.7_8+6.0_8)*(2_8*g+sqp) / (2.0_8*(g*g-1.0_8)+1e-16_8)
            val = (p1+p2)*cons
            x32 = tz*tz
            mp2 = SQRT(1.0_8-x32)
            sal_y = (sz*(1.0_8-x32)-tz*(tx*sx+ty*sy))/mp2*val
            sal_x = (tx*sy-ty*sx)/mp2*val
        END IF
    END SUBROUTINE

    SUBROUTINE proxy_source_compute(proxy_source_weights, tree_panels, points_panels, xs, ys, zs, sshs, areas, interp_degree, &
                                    point_count)
        ! compute the proxy source contributions from each point that the processor owns
        real(8), intent(out) :: proxy_source_weights(:)
        real(8), intent(in) :: xs(:), ys(:), zs(:), sshs(:), areas(:)
        type(cube_panel), intent(in) :: tree_panels(:)
        integer, intent(in) :: points_panels(:,:), interp_degree, point_count
        real(8) :: x, y, z, min_xi, max_xi, min_eta, max_eta, xi, eta, area, ssh
        integer :: i, j, i_t, k, l, shift, offset
        real(8), allocatable:: basis_vals(:,:)

        DO i = 1, point_count
            x = xs(i)
            y = ys(i)
            z = zs(i)
            area = areas(i)
            ssh = sshs(i)
            panelloop: DO j = 1, size(points_panels(:,i))
                ! loop over panels containing point i
                i_t = points_panels(j, i)
                IF (i_t == -1) THEN
                    exit panelloop
                ELSE
                    shift = (i_t-1)*(interp_degree+1)*(interp_degree+1)
                    min_xi = tree_panels(i_t)%min_xi
                    max_xi = tree_panels(i_t)%max_xi
                    min_eta = tree_panels(i_t)%min_eta
                    max_eta = tree_panels(i_t)%max_eta
                    call xieta_from_xyz(x, y, z, xi, eta, tree_panels(i_t)%face)
                    call interp_vals_bli(basis_vals, xi, eta, min_xi, max_xi, min_eta, max_eta, interp_degree)
                    offset = 0
                    DO k = 1, interp_degree+1 !eta loop
                        DO l = 1, interp_degree+1 ! xi loop
                            offset = offset+1
                            proxy_source_weights(shift+offset) = proxy_source_weights(shift+offset) + &
                                basis_vals(k, l)*area*ssh
                        END DO
                    END DO
                END IF
            END DO panelloop
        END DO
    END SUBROUTINE

    SUBROUTINE pc_interact(sal_longrad, sal_latgrad, i_t, source_panel, x, y, z, interp_degree, proxy_weights)
        real(8), intent(inout) :: sal_longrad(:), sal_latgrad(:)
        integer, intent(in) :: i_t, interp_degree
        real(8), intent(in) :: proxy_weights(:), x, y, z
        type(cube_panel), intent(in) :: source_panel
        real(8) :: min_xi, max_xi, min_eta, max_eta, eta, xi, cx, cy, cz, sal_x, sal_y
        real(8), allocatable :: cheb_xi(:), cheb_eta(:)
        integer :: i, j, offset

        min_xi = source_panel%min_xi
        max_xi = source_panel%max_xi
        min_eta = source_panel%min_eta
        max_eta = source_panel%max_eta

        ! barycentric lagrange interpolation points
        call bli_interp_points_shift(cheb_xi, min_xi, max_xi, interp_degree)
        call bli_interp_points_shift(cheb_eta, min_eta, max_eta, interp_degree)

        ! interact proxy source points and target point
        offset = 0
        DO j = 1, interp_degree+1
            eta = cheb_eta(j)
            DO i = 1, interp_degree+1
                xi = cheb_xi(i)
                call xyz_from_xieta(cx, cy, cz, xi, eta, source_panel%face)
                call sal_grad_gfunc(x, y, z, cx, cy, cz, sal_x, sal_y)
                offset = offset+1
                sal_longrad(i_t) = sal_longrad(i_t) + sal_x*proxy_weights(offset)
                sal_latgrad(i_t) = sal_latgrad(i_t) + sal_y*proxy_weights(offset)
            END DO
        END DO
    END SUBROUTINE

    SUBROUTINE pp_interact(sal_longrad, sal_latgrad, i_t, source_panel, xs, ys, zs, as, sshs, e_xs, e_ys, e_zs, e_as, &
                            e_sshs, own_points)
        real(8), intent(inout) :: sal_longrad(:), sal_latgrad(:)
        integer, intent(in) :: i_t, own_points
        real(8), intent(in) :: xs(:), ys(:), zs(:), as(:), sshs(:), e_xs(:), e_ys(:), e_zs(:), e_as(:), e_sshs(:)
        type(cube_panel), intent(in) :: source_panel
        real(8) :: tx, ty, tz, sx, sy, sz, ssh, area, sal_x, sal_y
        integer :: i, i_s, i_sr

        tx = xs(i_t)
        ty = ys(i_t)
        tz = zs(i_t)

        ! print *, i_t

        DO i = 1, source_panel%panel_point_count
            i_s = source_panel%relabeled_points_inside(i)
            ! i_sr = relabel_map(i_s)
            IF (i_s > own_points) THEN
                ! point is not own point
                i_s = i_s - own_points
                sx = e_xs(i_s)
                sy = e_ys(i_s)
                sz = e_zs(i_s)
                ssh = e_sshs(i_s)
                area = e_as(i_s)
            ELSE
                sx = xs(i_s)
                sy = ys(i_s)
                sz = zs(i_s)
                ssh = sshs(i_s)
                area = as(i_s)
            END IF
            call sal_grad_gfunc(tx, ty, tz, sx, sy, sz, sal_x, sal_y)
            ! print *, sal_x, sal_y
            ! print *, ssh, area
            ! print *, sal_latgrad(i_t)
            sal_longrad(i_t) = sal_longrad(i_t) + ssh*area*sal_x
            sal_latgrad(i_t) = sal_latgrad(i_t) + ssh*area*sal_y
            ! print *, sal_latgrad(i_t)
        END DO
    END SUBROUTINE

    ! computes all the particle-cluster interactions
    SUBROUTINE pc_compute(sal_longrad, sal_latgrad, pc_interaction_list, cube_panels, xs, ys, zs, interp_degree, &
                            proxy_weights)
        real(8), intent(in) :: xs(:), ys(:), zs(:), proxy_weights(:)
        type(interaction_pair), intent(in) :: pc_interaction_list(:)
        real(8), intent(inout) :: sal_longrad(:), sal_latgrad(:)
        type(cube_panel), intent(in) :: cube_panels(:)
        integer, intent(in) :: interp_degree
        integer :: i, i_s, proxy_count, i_t
        real(8), allocatable :: source_proxy_weights(:)

        proxy_count = (interp_degree+1)*(interp_degree+1)

        allocate(source_proxy_weights(proxy_count))

        DO i = 1, size(pc_interaction_list)
            i_s = pc_interaction_list(i)%index_source
            i_t = pc_interaction_list(i)%index_target
            source_proxy_weights = proxy_weights(proxy_count*(i_s-1)+1:proxy_count*i_s)
            call pc_interact(sal_longrad, sal_latgrad, i_t, cube_panels(i_s), xs(i_t), ys(i_t), zs(i_t), interp_degree, &
                                source_proxy_weights)
        END DO
    END SUBROUTINE

    ! computes all the particle particle interactions
    SUBROUTINE pp_compute(sal_longrad, sal_latgrad, pp_interaction_list, cube_panels, xs, ys, zs, areas, sshs, e_xs, &
                            e_ys, e_zs, e_areas, e_sshs, own_points)
        real(8), intent(inout) :: sal_longrad(:), sal_latgrad(:)
        type(interaction_pair), intent(in) :: pp_interaction_list(:)
        type(cube_panel), intent(in) :: cube_panels(:)
        real(8), intent(in) :: xs(:), ys(:), zs(:), areas(:), sshs(:), e_xs(:), e_ys(:), e_zs(:), e_areas(:), e_sshs(:)
        integer, intent(in) :: own_points
        integer :: i, i_t, i_s

        DO i = 1, size(pp_interaction_list)
            i_t = pp_interaction_list(i)%index_target
            i_s = pp_interaction_list(i)%index_source
            call pp_interact(sal_longrad, sal_latgrad, i_t, cube_panels(i_s), xs, ys, zs, areas, sshs, e_xs, e_ys, e_zs, &
                                e_areas, e_sshs, own_points)
        END DO
    END SUBROUTINE
END MODULE
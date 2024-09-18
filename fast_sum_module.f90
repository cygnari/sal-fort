MODULE FAST_SUM_MODULE
    use tree_module
    use bli_module
    use cube_sphere_module
    implicit none

    CONTAINS

    SUBROUTINE sal_grad_gfunc(tx, ty, tz, sx, sy, sz, sal_x, sal_y)
        real(8), intent(in) :: tx, ty, tz, sx, sy, sz
        real(8), intent(out) :: sal_x, sal_y
        real(8) :: g, mp, sqrtp, pi, cons, sqp, p1, p2, x32, val, mp2

        pi = 4.D0*DATAN(1.D0)
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
        real(8), allocatable :: cheb_xi(:), cheb_eta(:), basis_vals(:,:), proxy_weights(:,:)
        integer :: target_count, target_i, i, j, source_count, source_j, k
        real(8) :: sx, sy, sz, xi, eta, min_xi, max_xi, min_eta, max_eta, tx, ty, tz, cx, cy, cz, sal_x, sal_y

        target_count = target_panel%panel_point_count
        source_count = source_panel%panel_point_count

        min_xi = source_panel%min_xi
        max_xi = source_panel%max_xi
        min_eta = source_panel%min_eta
        max_eta = source_panel%max_eta

        ! barycentric lagrange interpolation points
        call bli_interp_points_shift(cheb_xi, min_xi, max_xi, interp_degree)
        call bli_interp_points_shift(cheb_eta, min_eta, max_eta, interp_degree)

        allocate(proxy_weights(interp_degree+1, interp_degree+1), source=0.0_8)

        ! compute proxy weights
        DO i = 1, source_count
            source_j = source_panel%points_inside(i)
            sx = xs(source_j)
            sy = ys(source_j)
            sz = zs(source_j)
            call xieta_from_xyz(sx, sy, sz, xi, eta, source_panel%face)
            call interp_vals_bli(basis_vals, xi, eta, min_xi, max_xi, min_eta, max_eta, interp_degree)
            DO k = 1, interp_degree+1 ! loop over eta
                DO j = 1, interp_degree+1 ! loop over xi
                    proxy_weights(j,k) = proxy_weights(j,k) + basis_vals(j, k)*area(source_j)*ssh(source_j)
                END DO
            END DO
        END DO

        ! interact proxy source points and target points
        DO i = 1, target_count
            target_i = target_panel%points_inside(i)
            tx = xs(target_i)
            ty = ys(target_i)
            tz = zs(target_i)
            DO k = 1, interp_degree+1
                eta = cheb_eta(k)
                DO j = 1, interp_degree+1
                    xi = cheb_xi(j)
                    call xyz_from_xieta(cx, cy, cz, xi, eta, source_panel%face)
                    call sal_grad_gfunc(tx, ty, tz, cx, cy, cz, sal_x, sal_y)
                    sal_longrad(target_i) = sal_longrad(target_i) + sal_x*proxy_weights(j,k)
                    sal_latgrad(target_i) = sal_latgrad(target_i) + sal_y*proxy_weights(j,k)
                END DO
            END DO
        END DO

    END SUBROUTINE

    SUBROUTINE cp_interact(sal_longrad, sal_latgrad, target_panel, source_panel, xs, ys, zs, area, ssh, interp_degree)
        real(8), intent(in) :: xs(:), ys(:), zs(:), area(:), ssh(:)
        type(cube_panel), intent(in) :: target_panel, source_panel
        integer, intent(in) :: interp_degree
        real(8), intent(inout) :: sal_longrad(:), sal_latgrad(:)
        real(8), allocatable :: cheb_xi(:), cheb_eta(:), sxs(:), sys(:), szs(:), areas(:), sshs(:), &
                                proxy_pots_x(:,:), proxy_pots_y(:,:), basis_vals(:,:)
        real(8) :: sx, sy, sz, tx, ty, tz, min_xi, max_xi, min_eta, max_eta, eta, cx, cy, cz, sal_x, sal_y, xi
        integer :: i, j, target_count, source_count, target_i, source_j, k

        target_count = target_panel%panel_point_count
        source_count = source_panel%panel_point_count

        min_xi = target_panel%min_xi
        max_xi = target_panel%max_xi
        min_eta = target_panel%min_eta
        max_eta = target_panel%max_eta

        call bli_interp_points_shift(cheb_xi, min_xi, max_xi, interp_degree)
        call bli_interp_points_shift(cheb_eta, min_eta, max_eta, interp_degree)

        allocate(proxy_pots_x(interp_degree+1, interp_degree+1), source=0.0_8)
        allocate(proxy_pots_y(interp_degree+1, interp_degree+1), source=0.0_8)

        allocate(sxs(source_count))
        allocate(sys(source_count))
        allocate(szs(source_count))
        allocate(areas(source_count))
        allocate(sshs(source_count))

        DO j = 1, source_count
            source_j = source_panel%points_inside(j)
            sxs(j) = xs(source_j)
            sys(j) = ys(source_j)
            szs(j) = zs(source_j)
            areas(j) = area(source_j)
            sshs(j) = ssh(source_j)
        END DO

        ! compute potentials at proxy targets
        DO j = 1, interp_degree+1 ! eta loop
            eta = cheb_eta(j)
            DO i = 1, interp_degree+1 ! xi loop
                call xyz_from_xieta(cx, cy, cz, cheb_xi(i), eta, target_panel%face)
                DO k = 1, source_count
                    call sal_grad_gfunc(cx, cy, cz, sxs(k), sys(k), szs(k), sal_x, sal_y)
                    proxy_pots_x(i,j) = proxy_pots_x(i,j) + sal_x*areas(k)*sshs(k)
                    proxy_pots_y(i,j) = proxy_pots_y(i,j) + sal_y*areas(k)*sshs(k)
                END DO
            END DO
        END DO

        ! interpolate from proxy targets to real targets
        DO i = 1, target_count
            target_i = target_panel%points_inside(i)
            tx = xs(target_i)
            ty = ys(target_i)
            tz = zs(target_i)
            call xieta_from_xyz(tx, ty, tz, xi, eta, target_panel%face)
            call interp_vals_bli(basis_vals, xi, eta, min_xi, max_xi, min_eta, max_eta, interp_degree)
            DO k = 1, interp_degree+1 ! eta loop
                DO j = 1, interp_degree+1 ! xi loop
                    sal_longrad(target_i) = sal_longrad(target_i) + basis_vals(j, k)*proxy_pots_x(j, k)
                    sal_latgrad(target_i) = sal_latgrad(target_i) + basis_vals(j, k)*proxy_pots_y(j, k)
                END DO
            END DO
        END DO

    END SUBROUTINE

    SUBROUTINE cc_interact(sal_longrad, sal_latgrad, target_panel, source_panel, xs, ys, zs, area, ssh, interp_degree)
        real(8), intent(in) :: xs(:), ys(:), zs(:), area(:), ssh(:)
        type(cube_panel), intent(in) :: target_panel, source_panel
        integer, intent(in) :: interp_degree
        real(8), intent(inout) :: sal_longrad(:), sal_latgrad(:)
        real(8), allocatable :: cheb_xi_s(:), cheb_eta_s(:), cheb_xi_t(:), cheb_eta_t(:)
        real(8), allocatable :: proxy_weights(:,:), proxy_pots_x(:,:), proxy_pots_y(:,:), basis_vals(:,:)
        integer :: target_count, source_count, i, j, k, l, source_j, target_i
        real(8) :: min_xi, max_xi, min_eta, max_eta, xi_t, eta_t, xi_s, eta_s, cxt, cyt, czt, cxs, cys, czs
        real(8) :: sal_x, sal_y, sx, sy, sz, tx, ty, tz

        target_count = target_panel%panel_point_count
        source_count = source_panel%panel_point_count

        min_xi = source_panel%min_xi
        max_xi = source_panel%max_xi
        min_eta = source_panel%min_eta
        max_eta = source_panel%max_eta

        call bli_interp_points_shift(cheb_xi_s, min_xi, max_xi, interp_degree)
        call bli_interp_points_shift(cheb_eta_s, min_eta, max_eta, interp_degree)

        allocate(proxy_weights(interp_degree+1, interp_degree+1), source=0.0_8)

        ! compute proxy weights
        DO i = 1, source_count
            source_j = source_panel%points_inside(i)
            sx = xs(source_j)
            sy = ys(source_j)
            sz = zs(source_j)
            call xieta_from_xyz(sx, sy, sz, xi_s, eta_s, source_panel%face)
            call interp_vals_bli(basis_vals, xi_s, eta_s, min_xi, max_xi, min_eta, max_eta, interp_degree)
            DO k = 1, interp_degree+1 ! loop over eta
                DO j = 1, interp_degree+1 ! loop over xi
                    proxy_weights(j,k) = proxy_weights(j,k) + basis_vals(j, k)*area(source_j)*ssh(source_j)
                END DO
            END DO
        END DO

        min_xi = target_panel%min_xi
        max_xi = target_panel%max_xi
        min_eta = target_panel%min_eta
        max_eta = target_panel%max_eta

        call bli_interp_points_shift(cheb_xi_t, min_xi, max_xi, interp_degree)
        call bli_interp_points_shift(cheb_eta_t, min_eta, max_eta, interp_degree)

        allocate(proxy_pots_x(interp_degree+1, interp_degree+1), source=0.0_8)
        allocate(proxy_pots_y(interp_degree+1, interp_degree+1), source=0.0_8)

        ! compute proxy potentials from proxy weights
        DO j = 1, interp_degree+1 ! target eta loop
            eta_t = cheb_eta_t(j)
            DO i = 1, interp_degree+1 ! target xi loop
                xi_t = cheb_xi_t(i)
                call xyz_from_xieta(cxt, cyt, czt, xi_t, eta_t, target_panel%face)
                DO l = 1, interp_degree+1 ! source eta loop
                    eta_s = cheb_eta_s(l)
                    DO k = 1, interp_degree+1 ! source xi loop
                        xi_s = cheb_xi_s(k)
                        call xyz_from_xieta(cxs, cys, czs, xi_s, eta_s, source_panel%face)
                        call sal_grad_gfunc(cxt, cyt, czt, cxs, cys, czs, sal_x, sal_y)
                        proxy_pots_x(i, j) = proxy_pots_x(i, j) + sal_x*proxy_weights(k, l)
                        proxy_pots_y(i, j) = proxy_pots_y(i, j) + sal_y*proxy_weights(k, l)
                    END DO
                END DO
            END DO
        END DO

        ! interpolate from proxy target points to real targets
        DO i = 1, target_count
            target_i = target_panel%points_inside(i)
            tx = xs(target_i)
            ty = ys(target_i)
            tz = zs(target_i)
            call xieta_from_xyz(tx, ty, tz, xi_t, eta_t, target_panel%face)
            call interp_vals_bli(basis_vals, xi_t, eta_t, min_xi, max_xi, min_eta, max_eta, interp_degree)
            DO k = 1, interp_degree+1 ! eta loop
                DO j = 1, interp_degree+1 ! xi loop
                    sal_longrad(target_i) = sal_longrad(target_i) + basis_vals(j, k)*proxy_pots_x(j, k)
                    sal_latgrad(target_i) = sal_latgrad(target_i) + basis_vals(j, k)*proxy_pots_y(j, k)
                END DO
            END DO
        END DO
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
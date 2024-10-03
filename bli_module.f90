MODULE BLI_MODULE
    implicit none

    CONTAINS

    SUBROUTINE bli_coeffs(coeffs, degree)
        real(8), allocatable, intent(out) :: coeffs(:)
        integer, intent(in) :: degree
        integer :: i

        allocate(coeffs(degree+1))

        IF (degree == 0) THEN
            coeffs(1) = 1.0_8
        ELSE
            DO i = 1, degree+1
                IF (i == 1) THEN
                    coeffs(i) = 0.5_8
                ELSE IF (i == degree + 1) THEN
                    coeffs(i) = 0.5_8 * ((-1) ** (i+1))
                ELSE
                    coeffs(i) = (-1.0_8) ** (i+1)
                END IF
            END DO
        END IF
    END SUBROUTINE

    SUBROUTINE bli_interp_points_shift(interp_points, min_x, max_x, degree)
        real(8), allocatable, intent(out) :: interp_points(:)
        real(8), intent(in) :: min_x, max_x
        integer, intent(in) :: degree
        real(8) :: x_range, x_shift, pi
        integer :: i

        pi = 4.D0*DATAN(1.D0)

        allocate(interp_points(degree+1))
        x_range = 0.5*(max_x - min_x)
        x_shift = 0.5*(max_x + min_x)

        IF (degree == 0) THEN
            interp_points(1) = x_shift
        ELSE 
            DO i = 1, degree + 1
                interp_points(i) = COS(pi/degree*(i-1))*x_range + x_shift
            END DO
        END IF
    END SUBROUTINE

    SUBROUTINE interp_vals_bli(vals, xi, eta, min_xi, max_xi, min_eta, max_eta, degree)
        real(8), allocatable, intent(out) :: vals(:,:)
        real(8), intent(in) :: xi, eta, min_xi, max_xi, min_eta, max_eta
        integer, intent(in) :: degree
        real(8), allocatable :: bli_xi_vals(:), bli_eta_vals(:), bli_coeff_vals(:), xi_func_vals(:), eta_func_vals(:)
        logical :: found_xi, found_eta
        integer :: i, j
        real(8) :: denom_xi, denom_eta, val

        allocate(vals(degree+1, degree+1))
        allocate(xi_func_vals(degree+1), source=0.0_8)
        allocate(eta_func_vals(degree+1), source=0.0_8)

        call bli_interp_points_shift(bli_xi_vals, min_xi, max_xi, degree)
        call bli_interp_points_shift(bli_eta_vals, min_eta, max_eta, degree)
        call bli_coeffs(bli_coeff_vals, degree)

        found_xi = .false.
        found_eta = .false.

        ! first check if xi/eta are an interpolation point
        DO i = 1, degree+1
            IF (ABS(xi - bli_xi_vals(i)) < 1.0e-16_8) THEN
                found_xi = .true.
                xi_func_vals(i) = 1.0_8
            END IF
            IF (ABS(eta - bli_eta_vals(i)) < 1.0e-16_8) THEN
                found_eta = .true.
                eta_func_vals(i) = 1.0_8
            END IF
        END DO

        ! xi/eta are not interpolation point, compute all the BLI basis values
        IF (.not. found_xi) THEN
            denom_xi = 0.0_8
            DO i = 1, degree+1
                val = bli_coeff_vals(i) / (xi - bli_xi_vals(i))
                xi_func_vals(i) = val
                denom_xi = denom_xi + val
            END DO
            DO i = 1, degree + 1
                xi_func_vals(i) = xi_func_vals(i) / denom_xi
            END DO
        END IF

        IF (.not. found_eta) THEN
            denom_eta = 0.0_8
            DO i = 1, degree+1
                val = bli_coeff_vals(i) / (eta - bli_eta_vals(i))
                eta_func_vals(i) = val
                denom_eta = denom_eta + val
            END DO
            DO i = 1, degree + 1
                eta_func_vals(i) = eta_func_vals(i) / denom_eta
            END DO
        END IF

        ! compute the outer product
        DO j = 1, degree + 1
            DO i = 1, degree + 1
                vals(i, j) = xi_func_vals(i) * eta_func_vals(j)
            END DO
        END DO
    END SUBROUTINE

END MODULE
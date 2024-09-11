MODULE CUBE_SPHERE_MODULE
    use mpi
    implicit none

    CONTAINS

    ! takes (x, y, z) coordinates of a point and returns which face of the cubed sphere it corresponds to
    INTEGER FUNCTION face_from_xyz(x, y, z) 
        real(8), intent(in) :: x, y, z
        real(8) :: ax, ay, az
        integer :: face
        ax = abs(x)
        ay = abs(y)
        az = abs(z)
        IF ((ax >= ay) .and. (ax >= az)) THEN
            IF (x >= 0) THEN
            face = 1
            ELSE 
            face = 3
            END IF
        ELSE IF ((ay >= ax) .and. (ay >= az)) THEN
            IF (y >= 0) THEN
            face = 2
            ELSE
            face = 4
            END IF
        ELSE
            IF (z >= 0) THEN
            face = 5
            ELSE
            face = 6
            END IF
        END IF
        face_from_xyz = face
        RETURN
    END FUNCTION

    ! computes xi eta angle coordinates on the cubed sphere from x, y, z points and optional face
    SUBROUTINE xieta_from_xyz(x, y, z, xi, eta, in_face)
        real(8), intent(in) :: x, y, z
        integer, intent(in), optional :: in_face
        real(8), intent(out) :: xi, eta
        integer :: face
        IF (.not. present(in_face)) THEN
            face = face_from_xyz(x, y, z)
        ELSE 
            face = in_face
        END IF
        IF (face == 1) THEN
            xi = ATAN(y/x); eta = ATAN(z/x)
        ELSE IF (face == 2) THEN
            xi = ATAN(-x/y); eta = ATAN(z/y)
        ELSE IF (face == 3) THEN
            xi = ATAN(y/x); eta = ATAN(-z/x)
        ELSE IF (face == 4) THEN
            xi = ATAN(-x/y); eta = ATAN(-z/y)
        ELSE IF (face == 5) THEN
            xi = ATAN(y/z); eta = ATAN(-x/z)
        ELSE IF (face == 6) THEN
            xi = ATAN(-y/z); eta = ATAN(-x/z)
        END IF
    END SUBROUTINE

    ! converts xi eta on the cubed sphere to x y z
    SUBROUTINE xyz_from_xieta(x, y, z, xi, eta, face) 
        real(8), intent(in) :: xi, eta
        integer, intent(in) :: face
        real(8), intent(out) :: x, y, z
        real(8) :: ax, ay
        ax = TAN(xi); ay = TAN(eta)
        IF (face == 1) THEN
            x = 1.0_8 / SQRT(1+ax*ax+ay*ay)
            y = ax*x
            z = ay*x
        ELSE IF (face == 2) THEN
            y = 1.0_8 / SQRT(1+ax*ax+ay*ay)
            x = -ax*y
            z = ay*y
        ELSE IF (face == 3) THEN
            x = -1.0_8 / SQRT(1+ax*ax+ay*ay)
            y = ax*x
            z = -ay*x
        ELSE IF (face == 4) THEN
            y = -1.0_8 / SQRT(1+ax*ax+ay*ay)
            x = -ax*y
            z = -ay*y
        ELSE IF (face == 5) THEN
            z = 1.0_8 / SQRT(1+ax*ax+ay*ay)
            x = -ay*z
            y = ax*z
        ELSE 
            z = -1.0_8 / SQRT(1+ax*ax+ay*ay)
            x = -ay*z
            y = -ax*z
        END IF
    END SUBROUTINE
END MODULE
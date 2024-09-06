program main
    use mpi
    implicit none

    real(8) theta, area, t1, t2
    integer point_count, ierr, process_rank, cluster_size, i
    real(8), allocatable :: areas(:), xs(:), ys(:), zs(:), sshs(:)
    NAMELIST /params/ theta, point_count

    call MPI_INIT(ierr)
    call MPI_COMM_SIZE(MPI_COMM_WORLD, cluster_size, ierr)
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

    IF (process_rank == 0) THEN
        open (file='./data-files/2313486_mpas_grid_areas.csv', unit=10) ! read in values from data files
        open (file='./data-files/2313486_mpas_grid_x.csv', unit=11)
        open (file='./data-files/2313486_mpas_grid_y.csv', unit=12)
        open (file='./data-files/2313486_mpas_grid_z.csv', unit=13)
        open (file='./data-files/2313486_mpas_grid_ssh.csv', unit=14)

        DO i = 1, point_count
            read (10, *) areas(i)
            read (11, *) xs(i)
            read (12, *) ys(i)
            read (13, *) zs(i)
            read (14, *) sshs(i)
        END DO

        close (10)
        close (11)
        close (12)
        close (13)
        close (14)
    END IF

    call Distribute_Array(areas, point_count, 0, MPI_COMM_WORLD, process_rank)
    call Distribute_Array(xs, point_count, 0, MPI_COMM_WORLD, process_rank)
    call Distribute_Array(ys, point_count, 0, MPI_COMM_WORLD, process_rank)
    call Distribute_Array(zs, point_count, 0, MPI_COMM_WORLD, process_rank)
    call Distribute_Array(sshs, point_count, 0, MPI_COMM_WORLD, process_rank)

    IF (process_rank == 0) THEN
        CALL CPU_TIME(t2)
        print *, 'run time : ',(t2 - t1)
    END IF

    IF (process_rank == 1) THEN
        DO i=1,point_count
            area = area + areas(i)
        END DO
        print *, "processor ", process_rank, " total area ", area
    END IF

    call MPI_FINALIZE(ierr)

end program

! Array A is fully located at processor source, distribute so that each processor has a copy of A
SUBROUTINE Distribute_Array(A, size, source, comm, p)
    use mpi
    integer size, source, comm, p, win, ierr
    real(8) :: A(size)
    integer (kind=MPI_ADDRESS_KIND) :: lowerbound, sizeofreal, disp
    disp = 0

    CALL MPI_TYPE_GET_EXTENT(MPI_DOUBLE_PRECISION, lowerbound, sizeofreal, ierr)
    CALL MPI_WIN_Create(A, size*sizeofreal, MPI_DOUBLE_PRECISION, MPI_INFO_NULL, comm, win, ierr)
    CALL MPI_WIN_Fence(0, win, ierr)
    IF (p /= source) THEN
        CALL MPI_Get(A, size, MPI_DOUBLE_PRECISION, 0, disp, size, MPI_DOUBLE_PRECISION, win, ierr)
    END IF
    CALL MPI_WIN_Fence(0, win, ierr)
    CALL MPI_WIN_Free(win, ierr)
END SUBROUTINE

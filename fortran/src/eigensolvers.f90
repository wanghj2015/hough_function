module eigensolvers

! Common interface to the symmetric eigenvalue solvers used to build Hough
! modes. F1/F2 (the matrices actually diagonalized here) are symmetric
! tridiagonal, so several standard library routines apply directly.
!
! All solve_symmetric() backends return eigenvalues/eigenvectors sorted in
! *descending* eigenvalue order, matching the convention of the Jacobi
! rotation method (Burkardt 2013) that the paper (Wang, Boyd & Akmaev 2016,
! Sect. 2.1) cites, and which this code uses as the default solver. Jacobi
! rotations and library routines (LAPACK's QR-based dsyev/dsyevd or the
! tridiagonal-specific dstev) agree on eigenvalues but pick different,
! equally valid sign gauges for the eigenvectors -- see docs/README.md.

implicit none
private
public :: solve_symmetric, compare_solvers, lapack_available, &
          solver_id_from_name, solver_name_from_id, &
          SOLVER_JACOBI, SOLVER_DSTEV, SOLVER_DSYEV, SOLVER_DSYEVD

integer, parameter :: SOLVER_JACOBI = 1
integer, parameter :: SOLVER_DSTEV  = 2
integer, parameter :: SOLVER_DSYEV  = 3
integer, parameter :: SOLVER_DSYEVD = 4

contains

  logical function lapack_available()
#ifdef HAVE_LAPACK
    lapack_available = .true.
#else
    lapack_available = .false.
#endif
  end function lapack_available


  integer function solver_id_from_name(name) result(id)
    character(len=*), intent(in) :: name

    select case (trim(adjustl(name)))
    case ('jacobi')
      id = SOLVER_JACOBI
    case ('dstev')
      id = SOLVER_DSTEV
    case ('dsyev')
      id = SOLVER_DSYEV
    case ('dsyevd')
      id = SOLVER_DSYEVD
    case default
      id = -1
    end select
  end function solver_id_from_name


  function solver_name_from_id(id) result(name)
    integer, intent(in) :: id
    character(len=16) :: name

    select case (id)
    case (SOLVER_JACOBI); name = 'jacobi'
    case (SOLVER_DSTEV);  name = 'dstev'
    case (SOLVER_DSYEV);  name = 'dsyev'
    case (SOLVER_DSYEVD); name = 'dsyevd'
    case default;         name = 'unknown'
    end select
  end function solver_name_from_id


  ! Solve A v = d v for a real symmetric (in fact tridiagonal) matrix A,
  ! returning eigenpairs sorted by descending eigenvalue.
  subroutine solve_symmetric(method, n, a, v, d)
    integer, intent(in)  :: method, n
    real(8), intent(in)  :: a(n,n)
    real(8), intent(out) :: v(n,n), d(n)

    real(8), allocatable :: awork(:,:), work(:), diag(:), offdiag(:)
    integer :: it_max, it_num, rot_num, info, lwork, i

    select case (method)
    case (SOLVER_JACOBI)
      allocate(awork(n,n))
      awork = a
      it_max = 1000
      call jacobi_eigenvalue(n, awork, it_max, v, d, it_num, rot_num)
      deallocate(awork)

    case (SOLVER_DSTEV)
#ifdef HAVE_LAPACK
      allocate(diag(n), offdiag(max(n-1,1)), work(max(1,2*n-2)))
      do i = 1, n
        diag(i) = a(i,i)
      end do
      do i = 1, n-1
        offdiag(i) = a(i,i+1)
      end do
      call dstev('V', n, diag, offdiag, v, n, work, info)
      call check_info('dstev', info)
      d = diag
      deallocate(diag, offdiag, work)
      call reverse_order(n, v, d)
#else
      call lapack_missing('dstev')
#endif

    case (SOLVER_DSYEV)
#ifdef HAVE_LAPACK
      v = a
      allocate(work(1))
      call dsyev('V', 'U', n, v, n, d, work, -1, info)
      lwork = int(work(1))
      deallocate(work)
      allocate(work(lwork))
      call dsyev('V', 'U', n, v, n, d, work, lwork, info)
      call check_info('dsyev', info)
      deallocate(work)
      call reverse_order(n, v, d)
#else
      call lapack_missing('dsyev')
#endif

    case (SOLVER_DSYEVD)
#ifdef HAVE_LAPACK
      v = a
      call solve_dsyevd(n, v, d)
      call reverse_order(n, v, d)
#else
      call lapack_missing('dsyevd')
#endif

    case default
      write(*,'(a,i0)') 'ERROR: unknown eigensolver id ', method
      stop 1
    end select

  end subroutine solve_symmetric


#ifdef HAVE_LAPACK
  subroutine solve_dsyevd(n, v, d)
    integer, intent(in)    :: n
    real(8), intent(inout) :: v(n,n)
    real(8), intent(out)   :: d(n)

    real(8), allocatable :: work(:)
    integer, allocatable :: iwork(:)
    integer :: lwork, liwork, info

    allocate(work(1), iwork(1))
    call dsyevd('V', 'U', n, v, n, d, work, -1, iwork, -1, info)
    lwork  = int(work(1))
    liwork = iwork(1)
    deallocate(work, iwork)

    allocate(work(lwork), iwork(liwork))
    call dsyevd('V', 'U', n, v, n, d, work, lwork, iwork, liwork, info)
    call check_info('dsyevd', info)
    deallocate(work, iwork)
  end subroutine solve_dsyevd
#endif


  subroutine reverse_order(n, v, d)
    ! LAPACK returns eigenvalues ascending; flip to descending so all
    ! solvers share the Jacobi routine's ordering convention.
    integer, intent(in)    :: n
    real(8), intent(inout) :: v(n,n), d(n)
    real(8) :: tmp(n)
    integer :: i

    d = d(n:1:-1)
    do i = 1, n/2
      tmp        = v(:,i)
      v(:,i)     = v(:,n-i+1)
      v(:,n-i+1) = tmp
    end do
  end subroutine reverse_order


  subroutine check_info(name, info)
    character(len=*), intent(in) :: name
    integer, intent(in) :: info
    if (info /= 0) then
      write(*,'(a,a,a,i0)') 'ERROR: ', name, ' failed with info = ', info
      stop 1
    end if
  end subroutine check_info


  subroutine lapack_missing(name)
    character(len=*), intent(in) :: name
    write(*,'(a,a,a)') 'ERROR: solver "', name, '" needs LAPACK, but this build has none.'
    write(*,'(a)')     '       Reconfigure with LAPACK available, or use --solver=jacobi.'
    stop 1
  end subroutine lapack_missing


  ! Solve the same matrix with every available solver and report how much
  ! their eigenvalues differ from the Jacobi reference (they should agree
  ! to machine precision -- eigenvalues are gauge-invariant even though
  ! eigenvector signs are not).
  subroutine compare_solvers(label, n, a)
    character(len=*), intent(in) :: label
    integer, intent(in) :: n
    real(8), intent(in) :: a(n,n)

    integer, parameter :: candidates(4) = &
        (/ SOLVER_JACOBI, SOLVER_DSTEV, SOLVER_DSYEV, SOLVER_DSYEVD /)
    real(8) :: v(n,n), d(n), d_ref(n)
    real(8) :: t0, t1
    integer :: i, m

    call solve_symmetric(SOLVER_JACOBI, n, a, v, d_ref)

    write(*,'(a)') ''
    write(*,'(a,a,a,i0,a)') ' Eigensolver comparison -- ', trim(label), ' (n = ', n, ')'
    write(*,'(a)') '   solver      max|lambda - lambda_jacobi|   time (s)'
    do i = 1, size(candidates)
      m = candidates(i)
      if (m /= SOLVER_JACOBI .and. .not. lapack_available()) cycle
      call cpu_time(t0)
      call solve_symmetric(m, n, a, v, d)
      call cpu_time(t1)
      write(*,'(3x,a10,4x,es16.6,4x,f8.4)') solver_name_from_id(m), maxval(abs(d - d_ref)), t1 - t0
    end do
    write(*,'(a)') '   (eigenvalues are gauge-invariant; eigenvector signs may differ'
    write(*,'(a)') '    between solvers -- Jacobi is used by default to match the published figures.)'

  end subroutine compare_solvers

end module eigensolvers

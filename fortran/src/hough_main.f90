program hough_main

! Compute Hough functions: scalar modes plus zonal/meridional wind modes,
! for a given zonal wavenumber s and normalized frequency sigma (= tidal
! frequency / (2*Omega)).
!
! Originally: October 2014 wanghoujun@gmail.com
! CLI / CMake / eigensolver-comparison rewrite: 2026

use eigensolvers, only: solver_id_from_name, solver_name_from_id, lapack_available, SOLVER_JACOBI

implicit none

real(8), parameter :: dpi = 3.1415926535897932d0, dtr = dpi/180.d0, rtd = 1.d0/dtr
real(8), parameter :: omega = 2*dpi/(24*3600.d0), earth_radius = 6.370d6
real(8), parameter :: g = 9.81d0

! ---- run configuration: defaults, overridable on the command line ----
integer            :: s = 2                ! zonal wavenumber
real(8)            :: sigma = 1.0d0        ! frequency / (2*Omega)
integer            :: nlat = 94            ! number of Gauss latitudes
integer            :: nn = 120             ! total modes (n2 = nn/2 per family)
integer            :: n2
character(len=16)  :: solver_name = 'jacobi'
integer            :: solver_id
character(len=16)  :: wind_method = 'auto'   ! auto | fd | groves
logical            :: use_fd
logical            :: compare_flag = .false.
logical            :: verbose = .false.
character(len=256) :: outfile = ''

real(8), allocatable :: latd(:), gx(:), gw(:)
real(8), allocatable :: theta(:,:), lambda(:)
real(8), allocatable :: theta_u(:,:), theta_v(:,:)
real(8), allocatable :: table1(:,:), table2(:,:)

integer :: i, j, k
real(8) :: wx
integer :: iounit

call parse_arguments()

n2 = nn/2
if (2*n2 /= nn) then
  write(*,'(a,i0,a)') 'ERROR: --nn must be even (got ', nn, ')'
  stop 1
end if

solver_id = solver_id_from_name(solver_name)
if (solver_id < 0) then
  write(*,'(a,a,a)') 'ERROR: unknown --solver "', trim(solver_name), '"'
  write(*,'(a)')     '       choose one of: jacobi, dstev, dsyev, dsyevd'
  stop 1
end if
if (solver_id /= SOLVER_JACOBI .and. .not. lapack_available()) then
  write(*,'(a,a,a)') 'ERROR: --solver=', trim(solver_name), ' needs LAPACK, but this build has none.'
  write(*,'(a)')     '       Reconfigure with LAPACK available, or use --solver=jacobi.'
  stop 1
end if

! Wind method: 'auto' picks fd for s=1 (where the Groves upward recurrence is
! numerically unstable) and groves otherwise; can be forced with --wind.
select case (trim(wind_method))
case ('auto')
  use_fd = (s == 1)
case ('fd')
  use_fd = .true.
case ('groves')
  use_fd = .false.
case default
  write(*,'(a,a,a)') 'ERROR: unknown --wind "', trim(wind_method), '" (choose auto, fd or groves)'
  stop 1
end select

if (len_trim(outfile) == 0) then
  write(outfile,'(a,i0,a,f0.4,a)') 'output/hough_s', s, '_sigma', sigma, '.dat'
end if

allocate(latd(nlat), gx(nlat), gw(nlat))
allocate(theta(nlat,nn), lambda(nn))
allocate(theta_u(nlat,nn), theta_v(nlat,nn))
allocate(table1(n2,n2), table2(n2,n2))

write(*,'(a)')           '=================================================================='
write(*,'(a)')           ' Hough function solver'
write(*,'(a,i0)')        '   zonal wavenumber s        = ', s
write(*,'(a,f0.4)')      '   normalized frequency sigma = ', sigma
write(*,'(a,i0)')        '   Gauss latitudes  (nlat)   = ', nlat
write(*,'(a,i0,a,i0,a)') '   modes per family (n2)     = ', n2, '  (', nn, ' total)'
write(*,'(a,a)')         '   eigensolver                = ', trim(solver_name)
if (use_fd) then
  write(*,'(a,a)')       '   wind method                = fd (differentiate scalar mode)'
else
  write(*,'(a,a)')       '   wind method                = groves (coefficient recurrence)'
end if
write(*,'(a,a)')         '   output file                = ', trim(outfile)
write(*,'(a)')           '=================================================================='

! Gauss-Legendre nodes and weights
call p_quadrature_rule(nlat, gx, gw)
latd = asin(gx)*rtd

if (verbose) then
  call r8vec2_print(nlat, gx, gw, '   Node      X              Weight')
  write(*,'(a)') ''
  write(*,'(a)') ' Gauss latitudes (deg):'
  write(*,'(8f9.3)') latd
end if

wx = sum(gw)
if (verbose) write(*,'(a,f12.8,a)') ' sum(weights) = ', wx, '  (expected 2)'

! Compute the scalar Hough modes
write(*,'(a)') ''
write(*,'(a)') ' -- scalar modes --'
call hough_mode(nlat, nn, n2, s, sigma, gx, gw, lambda, theta, solver_id, compare_flag)

! Compute the zonal/meridional wind Hough modes, either by differentiating
! the scalar modes (fd -- stable for all s) or via the Groves coefficient
! recurrence (groves -- fast, but unstable for s=1). The groves path
! re-solves the same eigenproblem internally, so the solver comparison table
! (if requested) is only printed once, above.
write(*,'(a)') ''
write(*,'(a)') ' -- wind (u,v) modes --'
if (use_fd) then
  call theta_to_theta_uv(nlat, nn, s, sigma, gx, theta, theta_u, theta_v)
else
  call hough_mode_uv(nlat, nn, n2, s, sigma, gx, gw, lambda, theta_u, theta_v, solver_id)
end if

! Sanity check: normalization of the u,v Hough modes
table1 = 0.d0
table2 = 0.d0
do i = 1, n2
  do j = 1, n2
    do k = 1, nlat
      table1(i,j) = table1(i,j) + gw(k)*theta_u(k,i)*theta_u(k,j)
      table2(i,j) = table2(i,j) + gw(k)*theta_v(k,i)*theta_v(k,j)
    end do
  end do
end do

write(*,'(a)') ''
write(*,'(a)') ' Normalization check (diagonal should be close to 1):'
write(*,'(a,es10.3)') '   max |norm(theta_u) - 1| = ', maxval(abs( (/ (table1(i,i), i=1,n2) /) - 1.d0 ))
write(*,'(a,es10.3)') '   max |norm(theta_v) - 1| = ', maxval(abs( (/ (table2(i,i), i=1,n2) /) - 1.d0 ))

write(*,'(a)') ''
write(*,'(a)') ' Gravest equivalent depths h = lambda * 4*(a*omega)^2/g  (km):'
write(*,'(a,5f10.4)') '   symmetric      : ', (lambda(i)*4*(earth_radius*omega)**2/g/1000.d0, i=1,min(5,n2))
write(*,'(a,5f10.4)') '   anti-symmetric : ', (lambda(n2+i)*4*(earth_radius*omega)**2/g/1000.d0, i=1,min(5,n2))

call execute_command_line('mkdir -p "' // trim(outdir(outfile)) // '"')

open(newunit=iounit, file=trim(outfile), form='unformatted', status='replace', action='write')
write(iounit) nlat, nn
write(iounit) latd
write(iounit) lambda
write(iounit) theta
write(iounit) theta_u
write(iounit) theta_v
close(iounit)

write(*,'(a)') ''
write(*,'(a,a)') ' wrote ', trim(outfile)
write(*,'(a)') '=================================================================='

contains

  function outdir(path) result(d)
    character(len=*), intent(in) :: path
    character(len=256) :: d
    integer :: p
    p = index(path, '/', back=.true.)
    if (p > 0) then
      d = path(1:p-1)
    else
      d = '.'
    end if
  end function outdir


  subroutine parse_arguments()
    integer :: nargs, iarg, eq
    character(len=256) :: arg, key, val

    nargs = command_argument_count()
    iarg = 1
    do while (iarg <= nargs)
      call get_command_argument(iarg, arg)

      select case (trim(arg))
      case ('-h', '--help')
        call print_help()
        stop
      case ('-v', '--verbose')
        verbose = .true.
      case ('--compare-solvers')
        compare_flag = .true.
      case default
        eq = index(arg, '=')
        if (eq > 0) then
          key = arg(1:eq-1)
          val = arg(eq+1:)
        else
          key = trim(arg)
          iarg = iarg + 1
          if (iarg > nargs) then
            write(*,'(a,a)') 'ERROR: missing value for ', trim(key)
            stop 1
          end if
          call get_command_argument(iarg, val)
        end if

        select case (trim(key))
        case ('-s', '--s')
          read(val,*) s
        case ('-f', '--sigma')
          read(val,*) sigma
        case ('--nlat')
          read(val,*) nlat
        case ('--nn')
          read(val,*) nn
        case ('--solver')
          solver_name = trim(val)
        case ('--wind')
          wind_method = trim(val)
        case ('-o', '--output')
          outfile = trim(val)
        case ('--preset')
          call apply_preset(trim(val))
        case default
          write(*,'(a,a,a)') 'ERROR: unrecognized option "', trim(key), '"'
          call print_help()
          stop 1
        end select
      end select

      iarg = iarg + 1
    end do
  end subroutine parse_arguments


  subroutine apply_preset(name)
    character(len=*), intent(in) :: name
    select case (name)
    case ('dw1')
      s = 1; sigma = 0.5d0
    case ('sw2')
      s = 2; sigma = 1.0d0
    case ('tw3')
      s = 3; sigma = 1.5d0
    case default
      write(*,'(a,a,a)') 'ERROR: unknown --preset "', name, '" (choose dw1, sw2 or tw3)'
      stop 1
    end select
  end subroutine apply_preset


  subroutine print_help()
    write(*,'(a)') 'Usage: hough_main [options]'
    write(*,'(a)') ''
    write(*,'(a)') 'Compute Hough functions (scalar + u/v wind modes) for zonal'
    write(*,'(a)') 'wavenumber s and normalized frequency sigma = frequency/(2*Omega).'
    write(*,'(a)') ''
    write(*,'(a)') 'Options:'
    write(*,'(a)') '  -s, --s=N            zonal wavenumber                    (default: 2)'
    write(*,'(a)') '  -f, --sigma=X        normalized frequency                 (default: 1.0)'
    write(*,'(a)') '      --preset=NAME    dw1 (s=1,sigma=0.5) | sw2 (s=2,sigma=1.0) | tw3 (s=3,sigma=1.5)'
    write(*,'(a)') '      --nlat=N         number of Gauss latitudes            (default: 94)'
    write(*,'(a)') '      --nn=N           total modes, must be even            (default: 120)'
    write(*,'(a)') '      --solver=NAME    jacobi | dstev | dsyev | dsyevd      (default: jacobi)'
    write(*,'(a)') '      --wind=METHOD    auto | fd | groves                   (default: auto)'
    write(*,'(a)') '                       auto = fd for s=1 (groves is unstable there), groves otherwise'
    write(*,'(a)') '      --compare-solvers  print an eigenvalue comparison across all available solvers'
    write(*,'(a)') '  -o, --output=PATH   output data file  (default: output/hough_s<S>_sigma<SIGMA>.dat)'
    write(*,'(a)') '  -v, --verbose        print extra diagnostic output'
    write(*,'(a)') '  -h, --help           show this help'
    write(*,'(a)') ''
    write(*,'(a)') 'Examples:'
    write(*,'(a)') '  hough_main --preset=dw1'
    write(*,'(a)') '  hough_main -s 2 -f 1.0 --solver=dstev --compare-solvers'
  end subroutine print_help

end program hough_main

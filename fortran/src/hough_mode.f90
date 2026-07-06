subroutine hough_mode(nlat, nn, n2, s, f, gx, gw, lambda, theta, solver, compare)

! Compute Hough modes
!
! November 2014 wanghoujun@gmail.com

implicit none

integer, intent(in) :: nlat, nn, n2

! eigensolver to use for F1/F2 (see eigensolvers.f90), and whether to print
! a comparison across all available solvers
integer, intent(in) :: solver
logical, intent(in) :: compare


real(8), parameter :: dpi = 3.1415926535, dtr = dpi/180., rtd = 1./dtr
real(8), parameter :: omega = 2*dpi/(24*3600.), a = 6.370e6
real(8), parameter :: g = 9.81, rd = 287., ts = 256., kappa = 2./7.


real(8), intent(in) :: f
integer, intent(in) :: s

! Expansion coefficients for normalized Hough functions and 
!                  normalized associated Legendre functions
real(8) :: v1(n2,n2), v2(n2,n2)

real(8) :: f1(n2,n2), f2(n2,n2)

real(8) :: d1(n2), d2(n2), lambda(nn)  ! eigenvalues

! Gauss-Legendre nodes and weights
real(8), intent(in), dimension(nlat) :: gx, gw

! Legendre polynomials
real(8) :: pnm(nlat, 0:nn+10)

! Hough modes: symmetric & anti-symmetric
real(8) :: theta1(nlat, n2), theta2(nlat, n2), theta(nlat, nn)

real(8) :: table1(n2,n2), table2(n2,n2)

integer :: i, j, k


! Hough expansion coefficients

call hough_coef(nn, s, f, f1, f2, d1, d2, v1, v2, solver, compare)

! Legendre functions

call pmn_polynomial ( nlat, nn+10, s, gx, pnm )

! Compute Hough modes

theta1 = 0.0d+00
theta2 = 0.0d+00

do i = 1, n2

do j = 1, n2
! Symmetric modes
if (s == 1) then
   k = 2*j-1
else if (s == 2) then
   k = 2*j
else if (s == 3) then
   k = 2*j+1
endif

theta1(:,i) = theta1(:,i) + v1(j,i)*pnm(:,k)

! Anti-symmetric modes
if (s == 1) then
   k = 2*j
else if (s == 2) then
   k = 2*j+1
else if (s == 3) then
   k = 2*j+2
endif

theta2(:,i) = theta2(:,i) + v2(j,i)*pnm(:,k)
enddo

enddo


! Check the normalization of the scalar Hough modes (should be ~1)

 table1 = 0.0d+00
 table2 = 0.0d+00
 do i = 1, n2
 do j = 1, n2
 do k = 1, nlat
 table1(i,j) = table1(i,j) + gw(k)*theta1(k,i)*theta1(k,j)
 table2(i,j) = table2(i,j) + gw(k)*theta2(k,i)*theta2(k,j)
 enddo
 enddo
 enddo

 write(*,'(a,es10.3)') '   max |norm(theta_sym)  - 1| = ', maxval(abs( (/ (table1(i,i), i=1,n2) /) - 1.d0 ))
 write(*,'(a,es10.3)') '   max |norm(theta_anti) - 1| = ', maxval(abs( (/ (table2(i,i), i=1,n2) /) - 1.d0 ))


! Put symmetric mode (theta1) and anti-symmetric mode (theta2) into theta
! and similarly for the associated eigenvalues d1 & d2 into lambda       

lambda(   1:n2) = d1(1:n2)
lambda(n2+1:nn) = d2(1:n2)

theta(:,   1:n2) = theta1(:,1:n2)
theta(:,n2+1:nn) = theta2(:,1:n2)


end subroutine hough_mode



subroutine hough_coef(n, s, f, f1, f2, d1, d2, v1, v2, solver, compare)

! Compute the coefficients of the Hough functions based on
! G. V. Groves (1981) Notes on obtaining the eigenvalues of
!                     Laplace tidal equation, Planet. Space Sci.
!
! August 2014 wanghoujun@gmail.com

use eigensolvers, only: solve_symmetric, compare_solvers
use, intrinsic :: ieee_arithmetic, only: ieee_is_finite

implicit none

integer, intent(in) :: n, s
!integer             :: n, s

real(8), intent(in) :: f
!real(8)             :: f

! matrix for symmetric & antisymmetric eigenvalues & eigenfunctions
real(8), intent(out) :: f1(n/2,n/2), f2(n/2,n/2)

! eigenvalues for symmetric & antisymmetric eigenfunctions
real(8), intent(out) :: d1(n/2), d2(n/2)
! symmetric & antisymmetric eigenfunctions
real(8), intent(out) :: v1(n/2,n/2), v2(n/2,n/2)

! eigensolver to use, and whether to print a solver comparison table
integer, intent(in) :: solver
logical, intent(in) :: compare

integer :: n2


real(8), parameter :: dpi = 3.1415926535, dtr = dpi/180., rtd = 1./dtr
real(8), parameter :: omega = 2*dpi/(24*3600.), a = 6.370e6
real(8), parameter :: g = 9.81, rd = 287., ts = 256., kappa = 2./7.


real(8) :: l(n), m(n)

real(8) :: sf, sigma
!real(8) :: sigma

!integer :: sf

integer :: i, j, r


! some constants

n2 = n/2

!f = 1.5d0
!f = 2.0d0
!f = 1.0d0
!f = 0.5d0

sigma = (2*f*omega)

sf = dble(dble(s)/f)

! define L(r) and M(r)
do r = s, n+s-1

i = r-s+1

! define L(r)
l(i) = dsqrt(dble((r+s+1)*(r+s+2)*(r-s+1)*(r-s+2)))/ &
       dble((2*r+3)*dsqrt(dble((2*r+1)*(2*r+5)))*dble(sf-(r+1)*(r+2)))


! define M(r)
if (s == 2 .and. r == 2) then
!m(i) = -dble(f**2-1)/dble((s/f+r)*(sf-r-1)) &
!       +dble((r-s+1)*(r+s+1)*(sf+r+2)) &
!       /dble((2*r+1)*(2*r+3)*(sf-r-1)*(sf-(r+1)*(r+2)))

 m(i) = -dble(f**2*(sf-r*(r+1)))/dble((r*(r+1))**2) &
        +dble((r+2)**2*(r+s+1)*(r-s+1)) &
        /dble((r+1)**2*(2*r+3)*(2*r+1)*(sf-(r+1)*(r+2))) 

!else if (s == 1 .and. r == 2) then
!m(i) = -dble(f**2*(sf-r*(r+1)))/dble((r*(r+1))**2) &
!       +dble((r+2)**2*(r+s+1)*(r-s+1)) &
!       /dble((r+1)**2*(2*r+3)*(2*r+1)*(sf-(r+1)*(r+2)))

!m(i) = -dble(f**2-1)/dble((sf+r)*(sf-r-1)) &
!       +dble((r-s+1)*(r+s+1)*(sf+r+2)) &
!       /dble((2*r+1)*(2*r+3)*(sf-r-1)*(sf-(r+1)*(r+2)))

!m(i) = -dble(f**2-1)/dble((sf+r)*(sf-r-1)) &
!       +dble((r-s)*(r+s)*(sf-r+1)) &
!       /dble((2*r-1)*(2*r+1)*(sf+r)*(sf+epsilon(sf)-r*(r-1))) &
!       +dble((r-s+1)*(r+s+1)*(sf+r+2)) &
!       /dble((2*r+1)*(2*r+3)*(sf-r-1)*(sf-(r+1)*(r+2)))

!elseif (s == 1 .and. r == 1) then
!m(i) = -dble(f**2*(sf-r*(r+1)))/dble((r*(r+1))**2) &
!       +dble((r+2)**2*(r+s+1)*(r-s+1)) &
!       /dble((r+1)**2*(2*r+3)*(2*r+1)*(sf-(r+1)*(r+2))) &
!       +dble((r-1)**2*(r**2-s**2)) &
!       /dble(r**2*(4*r**2-1)*(sf-r*(r-1)))

!else if (s == 2) then
!m(i) = -dble(f**2-1)/dble((sf+r)*(sf-r-1)) &
!       +dble((r-s)*(r+s)*(sf-r+1)) &
!       /dble((2*r-1)*(2*r+1)*(sf+r)*(sf-r*(r-1))) &
!       +dble((r-s+1)*(r+s+1)*(sf+r+2)) &
!       /dble((2*r+1)*(2*r+3)*(sf-r-1)*(sf-(r+1)*(r+2)))

!else if (s == 1) then
!m(i) = -dble(f**2*(sf-r*(r+1)))/dble((r*(r+1))**2) &
!       +dble((r+2)**2*(r+s+1)*(r-s+1)) &
!       /dble((r+1)**2*(2*r+3)*(2*r+1)*(sf-(r+1)*(r+2))) &
!       +dble((r-1)**2*(r**2-s**2)) &
!       /dble(r**2*(4*r**2-1)*(sf-r*(r-1)))

else

 m(i) = -dble(f**2*(sf-r*(r+1)))/dble((r*(r+1))**2) &
        +dble((r+2)**2*(r+s+1)*(r-s+1)) &
        /dble((r+1)**2*(2*r+3)*(2*r+1)*(sf-(r+1)*(r+2))) &
        +dble((r-1)**2*(r**2-s**2)) &
        /dble(r**2*(4*r**2-1)*(sf-r*(r-1)))


endif

! A few (s, f) combinations have a genuine infinite equivalent-depth mode
! at this r (an expected physical feature -- the paper's "[0]" missing
! mode -- not a bug; see docs/README.md). Clamp to the largest finite
! double, matching the Python port's identical guard, so every eigensolver
! can operate on a literal matrix; the eigenvalue/depth still overflows to
! Infinity downstream in the usual floating-point way.
if (.not. ieee_is_finite(m(i))) m(i) = huge(m(i))

enddo

! build F1 (odd r) and F2 (even r), the symmetric tridiagonal matrices
! whose eigenvalues/eigenvectors give the symmetric/anti-symmetric Hough modes
f1 = 0.
f2 = 0.

do i = 1, n2
f1(i,i) = m(2*i-1)
if (i+1 <= n2) then
   f1(i,i+1) = l(2*i-1)
   f1(i+1,i) = l(2*i-1)
endif
enddo

do i = 1, n2
f2(i,i) = m(2*i)
if (i+1 <= n2) then
   f2(i,i+1) = l(2*i)
   f2(i+1,i) = l(2*i)
endif
enddo


! F1, F2 are symmetric tridiagonal; solve with the requested eigensolver
! (default: Jacobi rotations, matching the paper's sign convention -- see
! docs/README.md). Optionally print a cross-solver comparison.

call solve_symmetric(solver, n2, f1, v1, d1)
if (compare) call compare_solvers('F1 (symmetric family)', n2, f1)

call solve_symmetric(solver, n2, f2, v2, d2)
if (compare) call compare_solvers('F2 (anti-symmetric family)', n2, f2)

end subroutine hough_coef



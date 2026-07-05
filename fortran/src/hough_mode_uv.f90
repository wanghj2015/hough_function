subroutine hough_mode_uv(nlat, nn, n2, s, f, gx, gw, lambda, theta_u, theta_v, solver)

! Compute Hough modes for wind
!
! February 2015 wanghoujun@gmail.com

implicit none

integer, intent(in) :: nlat, nn, n2

! eigensolver to use for F1/F2 (see eigensolvers.f90); the solver comparison
! table is not re-printed here since hough_mode() already solves and prints
! it for the same (s,f)
integer, intent(in) :: solver


real(8), parameter :: dpi = 3.1415926535, dtr = dpi/180., rtd = 1./dtr
real(8), parameter :: omega = 2*dpi/(24*3600.), a = 6.370e6
real(8), parameter :: g = 9.81, rd = 287., ts = 256., kappa = 2./7.


real(8), intent(in) :: f
integer, intent(in) :: s

! Expansion coefficients for normalized Hough functions and 
!                  normalized associated Legendre functions
real(8) :: a1(n2,n2), a2(n2,n2)

real(8) :: u1(n2,n2), u2(n2,n2)
real(8) :: v1(n2,n2), v2(n2,n2)

real(8) :: f1(n2,n2), f2(n2,n2)

real(8) :: d1(n2), d2(n2), lambda(nn)  ! eigenvalues

! Gauss-Legendre nodes and weights
real(8), intent(in), dimension(nlat) :: gx, gw

! Legendre polynomials
real(8) :: pnm(nlat, 0:nn+5)

! Hough modes: symmetric & anti-symmetric
real(8), dimension(nlat,n2) :: theta_u1, theta_u2, theta_v1, theta_v2
real(8), dimension(nlat,nn) :: theta_u, theta_v

real(8) :: table1(n2,n2), table2(n2,n2)

real(8) :: latd(nlat)

integer :: i, j, k


! Hough expansion coefficients (re-solves F1/F2; solver comparison is
! suppressed here since hough_mode() already printed it for this (s,f))

call hough_coef(nn, s, f, f1, f2, d1, d2, a1, a2, solver, .false.)

call hough_coef_uv(nn, s, f, a1, a2, u1, u2, v2, v1)

! Legendre functions

call pmn_polynomial ( nlat, nn+5, s, gx, pnm )

! Compute Hough modes

theta_u1 = 0.0d+00
theta_v1 = 0.0d+00
theta_u2 = 0.0d+00
theta_v2 = 0.0d+00

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

theta_u1(:,i) = theta_u1(:,i) + u1(j,i)*pnm(:,k)
theta_v1(:,i) = theta_v1(:,i) + v1(j,i)*pnm(:,k)

! Anti-symmetric modes
if (s == 1) then
   k = 2*j
else if (s == 2) then
   k = 2*j+1
else if (s == 3) then
   k = 2*j+2
endif

theta_u2(:,i) = theta_u2(:,i) + u2(j,i)*pnm(:,k)
theta_v2(:,i) = theta_v2(:,i) + v2(j,i)*pnm(:,k)
enddo

enddo


latd = asin(gx)*rtd

do k = 1, nlat
theta_u1(k,:) = theta_u1(k,:)/cos(latd(k)*dtr)
theta_u2(k,:) = theta_u2(k,:)/cos(latd(k)*dtr)
theta_v1(k,:) = theta_v1(k,:)/cos(latd(k)*dtr)
theta_v2(k,:) = theta_v2(k,:)/cos(latd(k)*dtr)
enddo


! Check the normalization of the u-wind Hough modes (should be ~1)

 table1 = 0.0d+00
 table2 = 0.0d+00
 do i = 1, n2
 do j = 1, n2
 do k = 1, nlat
 table1(i,j) = table1(i,j) + gw(k)*theta_u1(k,i)*theta_u1(k,j)
 table2(i,j) = table2(i,j) + gw(k)*theta_u2(k,i)*theta_u2(k,j)
 enddo
 enddo
 enddo

 write(*,'(a,es10.3)') '   max |norm(theta_u_sym)  - 1| = ', maxval(abs( (/ (table1(i,i), i=1,n2) /) - 1.d0 ))
 write(*,'(a,es10.3)') '   max |norm(theta_u_anti) - 1| = ', maxval(abs( (/ (table2(i,i), i=1,n2) /) - 1.d0 ))


! Put symmetric mode (theta1) and anti-symmetric mode (theta2) into theta
! and similarly for the associated eigenvalues d1 & d2 into lambda       

lambda(   1:n2) = d1(1:n2)
lambda(n2+1:nn) = d2(1:n2)

theta_u(:,   1:n2) = theta_u1(:,1:n2)
theta_u(:,n2+1:nn) = theta_u2(:,1:n2)

!theta_v(:,   1:n2) = theta_v1(:,1:n2)
!theta_v(:,n2+1:nn) = theta_v2(:,1:n2)

 theta_v(:,   1:n2) = theta_v2(:,1:n2)
 theta_v(:,n2+1:nn) = theta_v1(:,1:n2)


end subroutine hough_mode_uv



subroutine hough_coef_uv(n, s, f, a1, a2, u1, u2, v1, v2)

! Compute the coefficients of the Hough functions for wind based on 
! G. V. Groves (1981) Notes on obtaining the eigenvalues of
!                     Laplace tidal equation, Planet. Space Sci.
!
! February 2015 wanghoujun@gmail.com

implicit none

integer, intent(in) :: n, s
!integer             :: n, s

real(8), intent(in) :: f
!real(8)             :: f

! symmetric & antisymmetric eigenfunctions for temperature
real(8), dimension(n/2,n/2), intent(in) :: a1, a2

! symmetric & antisymmetric eigenfunctions for velocity u & v
real(8), dimension(n/2,n/2), intent(out) :: u1, u2, v1, v2

real(8), dimension(n,n) :: a, u, v

real(8), dimension(n) :: e
 
real(8), dimension(n/2) :: e1, e2


integer :: n2


real(8), parameter :: dpi = 3.1415926535, dtr = dpi/180., rtd = 1./dtr
!real(8), parameter :: omega = 2*dpi/(24*3600.), a = 6.370e6
real(8), parameter :: omega = 2*dpi/(24*3600.), re = 6.370e6
real(8), parameter :: g = 9.81, rd = 287., ts = 256., kappa = 2./7.


real(8) :: sf, sigma

!integer :: sf

integer :: i, j, r, i1, j1, i2, j2 


! some constants

n2 = n/2

!sigma = 2*f*omega
sigma = f

sf = dble(dble(s)/f)

! define e_r
do r = s, n+s-1
i = r-s+1
e(i) = dsqrt(dble(dble(r**2-s**2)/dble(4*r**2-1)))
enddo

do r = s, n+s-2, 2
i = ((r-s+1)+1)/2
e1(i) = dsqrt(dble(dble(r**2-s**2)/dble(4*r**2-1)))
enddo

do r = s+1, n+s-1, 2
i = (r-s+1)/2
e2(i) = dsqrt(dble(dble(r**2-s**2)/dble(4*r**2-1)))
enddo


! put the symmetric and anti-symmetric eigenfunctions into one matrix

a = 0.d0

do i = 1, n-1, 2
do j = 1, n-1, 2
i1 = (i+1)/2
j1 = (j+1)/2
a(i,j) = a1(i1,j1)
enddo
enddo

do i = 2, n, 2
do j = 2, n, 2
i2 = i/2
j2 = j/2
a(i,j) = a2(i2,j2)
enddo
enddo


! compute u_r & v_r recursively

u = 0.d0
v = 0.d0

u(1,:) = 1/sigma*(s+1)*e(1)*a(1,:)
v(1,:) = 0.d0

u(2,:) = (sigma*v(1,:) + 1/sigma*(s+2)*e(2)*a(2,:))/e(2)
v(2,:) = (sigma*u(1,:) - s/sigma*a(1,:))/e(2)

do r = s+2, n+s-1

i = r-s+1

u(i,:) = (sigma*v(i-1,:) - e(i-1)*u(i-2,:) &
       - 1/sigma*((r-2)*e(i-1)*a(i-2,:)-(r+1)*e(i)*a(i,:)))/e(i)

v(i,:) = (sigma*u(i-1,:) - e(i-1)*v(i-2,:) - s/sigma*a(i-1,:))/e(i)

enddo

! put the symmetric and anti-symmetric eigenfunctions into one matrix

u1 = 0.d0
u2 = 0.d0
v1 = 0.d0
v2 = 0.d0

do i = 1, n-1, 2
do j = 1, n-1, 2
i1 = (i+1)/2
j1 = (j+1)/2
u1(i1,j1) = u(i,j)
v2(i1,j1) = v(i,j)
enddo
enddo

do i = 2, n, 2
do j = 2, n, 2
i2 = i/2
j2 = j/2
u2(i2,j2) = u(i,j)
v1(i2,j2) = v(i,j)
enddo
enddo


! for symmetric modes

u1(1,:) = 1/sigma*(s+1)*a1(1,:) 
v1(1,:) = (sigma*u1(1,:) - s/sigma*a1(1,:))/e2(1)

!u1(2,:) = (sigma*v1(1,:) - e2(1)*u1(1,:) &
!         - 1/sigma*(s*e2(1)*a1(1,:) - (s+3)*e1(2)*a1(2,:)))/e1(2)
!v1(2,:) = (sigma*u1(2,:) - e1(2)*v1(1,:) - s/sigma*a1(2,:))/e2(2)

do r = s+2, n+s-1

if (mod(r-s,2) == 0) then

i = ((r-s+1)+1)/2
u1(i,:) = (sigma*v1(i-1,:) - e2(i-1)*u1(i-1,:) &
        - 1/sigma*((r-2)*e2(i-1)*a1(i-1,:)-(r+1)*e1(i)*a1(i,:)))/e1(i)

else

i = (r-s+1)/2
v1(i,:) = (sigma*u1(i,:) - e1(i)*v1(i-1,:) - s/sigma*a1(i,:))/e2(i)

endif

enddo


! for anti-symmetric modes

v2(1,:) = 0.d0
u2(1,:) = (sigma*v2(1,:) + 1/sigma*(s+2)*e2(1)*a2(1,:))/e2(1)

do r = s+2, n+s-1

if (mod(r-s,2) == 0) then

i = ((r-s+1)+1)/2
v2(i,:) = (sigma*u2(i-1,:) - e2(i-1)*v2(i-1,:) - s/sigma*a2(i-1,:))/e1(i)

else

i = (r-s+1)/2
u2(i,:) = (sigma*v2(i,:) - e1(i)*u2(i-1,:) &
        - 1/sigma*((r-2)*e1(i)*a2(i-1,:)-(r+1)*e2(i)*a2(i,:)))/e2(i)

endif

enddo


end subroutine hough_coef_uv


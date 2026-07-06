subroutine theta_to_theta_uv(nlat, nn, s, sigma, mu, theta, theta_u, theta_v)

! Compute the zonal/meridional wind Hough modes from the (already computed)
! scalar Hough modes by DIFFERENTIATING them -- the "fd" (finite-difference)
! method. This is an alternative to the Groves (1981) upward coefficient
! recurrence in hough_mode_uv(): the recurrence is numerically unstable for
! s=1 (it amplifies the high-degree eigenvector tail into a grid-scale
! sawtooth), whereas differentiating the well-resolved scalar mode is stable.
!
! Relations (same as the Python port's nalp_hough._winds; x = mu = sin(lat)):
!   b2 = sqrt(1-x^2)/(sigma^2 - x^2)
!   u  = b2 * ( s/(1-x^2) * theta        - x/sigma * d(theta)/dx )
!   v  = b2 * ( s/sigma * x/(1-x^2) * theta -        d(theta)/dx )
! The derivative is a non-uniform central difference on the Gauss grid.
!
! Ported from tidal_analysis/src/hough/theta_to_theta_uv.f90 (Houjun Wang,
! Feb 2015), corrected to use the computed derivative dxdl (the archived
! copy substituted b2 for the derivative in the second term).

implicit none

integer, intent(in) :: nlat, nn, s
real(8), intent(in) :: sigma
real(8), intent(in) :: mu(nlat)            ! = sin(latitude) = Gauss nodes gx
real(8), intent(in) :: theta(nlat,nn)      ! scalar Hough modes (cols = modes)
real(8), dimension(nlat,nn), intent(out) :: theta_u, theta_v

real(8) :: dxdl(nlat), b2
integer :: i, j, jp, jm

theta_u = 0.d0
theta_v = 0.d0

do i = 1, nn

  ! d(theta)/dmu via non-uniform central differences (one-sided at the ends)
  do j = 1, nlat
    jp = min(j+1, nlat)
    jm = max(j-1, 1)
    dxdl(j) = (theta(jp,i) - theta(jm,i)) / (mu(jp) - mu(jm))
  end do

  do j = 1, nlat
    b2 = sqrt(1.d0 - mu(j)**2) / (sigma**2 - mu(j)**2)
    theta_u(j,i) = b2 * ( dble(s)/(1.d0 - mu(j)**2)*theta(j,i) &
                          - mu(j)/sigma*dxdl(j) )
    theta_v(j,i) = b2 * ( dble(s)/sigma*mu(j)/(1.d0 - mu(j)**2)*theta(j,i) &
                          - dxdl(j) )
  end do

end do

end subroutine theta_to_theta_uv

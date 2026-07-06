cmake_minimum_required(VERSION 3.16)
project(hough_function Fortran)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

set(CMAKE_Fortran_FLAGS_RELEASE "-O2")
set(CMAKE_Fortran_FLAGS_DEBUG "-g -fbounds-check -Wall")

# .mod files clutter whatever directory the compiler runs in; keep them
# out of the source tree
set(CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/mod)

# LAPACK is optional: the Jacobi solver always works, but dstev/dsyev/
# dsyevd need LAPACK. On macOS this typically resolves via the Accelerate
# framework or Homebrew's openblas; on Linux, apt/yum's liblapack-dev.
find_package(LAPACK)

add_executable(hough_main
  src/hough_main.f90
  src/hough_mode.f90
  src/hough_mode_uv.f90
  src/theta_to_theta_uv.f90
  src/eigensolvers.f90
  src/jacobi_eigenvalue.f90
  src/legendre_polynomial_sub.f90
)

# -cpp turns on the C preprocessor so eigensolvers.f90 can #ifdef HAVE_LAPACK
target_compile_options(hough_main PRIVATE -cpp)

if(LAPACK_FOUND)
  message(STATUS "LAPACK found (${LAPACK_LIBRARIES}) -- dstev/dsyev/dsyevd solvers enabled")
  target_compile_definitions(hough_main PRIVATE HAVE_LAPACK)
  target_link_libraries(hough_main PRIVATE ${LAPACK_LIBRARIES})
else()
  message(STATUS "LAPACK not found -- only the Jacobi solver will be available (this is fine; it's the default)")
endif()

install(TARGETS hough_main RUNTIME DESTINATION bin)

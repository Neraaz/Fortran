!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! Let's diagonalize matrix A to get eigenvalues (l) and eigenvectors (x). Ax = lx
! (A-lI)x=0.
!
! Definitions
!
! A = transpose(A) ==> symmetric
! A = complex_conjugate(transpose(A)) ==> Hermitian (or Hermitian conjugate) or self-adjoint
! transpose(A).A = A.transpose(A) = I ==> orthogonal
! Hermitian_conjugate(A).A = A.Hermitian_conjugate(A) ==> Normal. Unitary if it
! equals to unitary matrix (I)
! 
! Matrix           Eigenvalues
! Real symmetric     real
! Hermitian          real
! Real nonsymmetric  may have both real and imaginary 
! non Hermitian 
! complex matrix     complex
!
!It looks like a difficult task. EISPACK is a package of fortran subroutine
!that perform diagonalization of matrix which is overtaken by LAPACK.
!
!LAPACK = LINPACK (solving linear equations AX=B, linear least square problems)
!+ EISPACK (solving matrix to get eigenvalues and eigenvectors)
!
!LAPACK routines are written so that as much as possible of the computation is
!performed by calls to the BLAS (Basic linear algebra subprograms)
!
! One can use LAPACK routine to pertorm those operations.
! Download lapack routine from http://www.netlib.org/lapack-3.10.0/
! 
! Extract it and install it using any compiler. We can do it using gnu compiler
! gcc/9.3.0. Compiling with gcc/4.8 creates *.a files, but doesn't execute the
! code.
!
! Copy make.inc.gfortran to make.inc (to compile using gfortran). 
! cp INSTALL/make.inc.gfortran make.inc (try with other compiler such as ifort,
! mpiifort, and mpifort.
!
!BLAS: Low level routine to perform linear algebra operations such as vector
!addition, scalar multiplication, dot products, linear combinations, and matrix
!multiplication.
!
!We can further customize according to our needs.
! Now simply type "make" inside lapack parent folder. This will create libraries in ".a" format.
! mv librefblas.a libblas.a
!
!Now we can link those libraries while compiling our simple fortran code.
!For dynamic linking: gfortran -o a.out code.f90 -L/path_to_our_lib.a_folder
!-llapack -lblas .... (other libraries if needed)
!
!For static linking: gfortran -o a.out code.f90 /path_to_lapack_lib/liblapack.a
!/path_to_blas_lib/libblas.a
!
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!Solving Linear equations AX=B
!
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
program main
implicit none
external :: dgesv
real(kind=8) :: a(2, 2) ! Matrix A.
real(kind=8) :: b(2) ! Vector b/x.
real(kind=8) :: pivot(2) ! Pivot indices (list of swap operations).
integer :: rc ! Return code.
a = reshape([ 2.d0, 3.d0, 1.d0, 1.d0 ], [ 2, 2 ])
b = [ 5.d0, 6.d0 ]
call dgesv(2, 1, a, 2, pivot, b, 2, rc)
if (rc /= 0) then
print '(a, i0)', 'Error: ', rc
stop
end if
print '("Solution (x1, x2): ", f0.4, ", ", f0.4)', b
end program main

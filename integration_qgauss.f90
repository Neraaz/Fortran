!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!Integration of a function
!Written by Niraj K. Nepal, adapted from Numerical Recipe: Fortran 90
!
!Usage: gfortran -o main integration_qgauss.f90 -L/path_to_lapack_libs -llapack -lblas
!./main
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

program main
implicit double precision (a-h,o-z)
integer i
real(kind=8) xmat(10,10),eig(10),wi(10),xvec(10,10)
real(kind=8) rslt
call GWmat(10,xmat,'gleg')
call diagonalize(10,100,xmat,eig,wi,xvec)
!call GolWel(10,xt,wt,'gleg')
!real(kind=8) rslt1, rslt2, rslt3
!a = aj(1,'gleg')
!external func
!call ser(1.d0,2.0d0,1000,f)
!write(*,*) 'array: ',f
do i=1,10   !printing grid points and corresponding weights.
write(*,*) "", eig(i),2.d0*xvec(1,i)**2.d0
enddo
call qgaus(1.d0,2.d0,10,rslt,'gleg')
write(*,*), "Integration: ", rslt
return
end program

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!function to integrate, Define your function here.
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
Double Precision Function func(x)
implicit double precision (a-h,o-z)
func=dsin(x)
return
end

!c ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
!Subroutine to perform integration.
!
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine qgaus(x1,x2,n,rslt,method)
implicit none
real(kind=8) x1,x2,rslt,func
CHARACTER(LEN=*), INTENT(IN) :: method
external func
integer j,n
real(kind=8) dx,xm,xr,w(n),x(n),xt,f
if (method .eq. 'gleg') then
 call GolWel(n,x,w,'gleg')
endif
xm=0.5d0*(x2+x1)
xr=0.5d0*(x2-x1)
rslt=0.d0
do j=1,n
  dx=xr*x(j)
  xt=xm+dx
  f=func(xt)
  rslt=rslt + xr*w(j)*f
enddo
return
end

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!Generating Gauss rule using Golub Welsch method. Any set of polynomials
!{Pj(x)}_(j=0,N-1) satisfies a 3 term recurrence relation. 
! Pj+1(x) = (aj*x + bj)Pj(x) - cj*Pj-1(x) with P0(x) = 1 and P1(x) = A0*x + B0
! Where aj, bj, and cj are defined as (https://dlmf.nist.gov/18.9)
! for j=0 to N-1
! Table
!Polynomial                 aj              bj              cj     mu0
!  type
!Gauss-legendre(gleg)     (2*j+1)/(j+1)    0             j/(j+1)   2.0
!
!Gauss-Laguerre(glag)     -1 / (j + 1)   (2j+alpha+1)/(j+1)  (j+alpha)/(j+1) 
!
!Let's create a symmetric tridiagonal matrix (J).
!  * * * * * * * * * * * * * * * * * * * * *
!  * al1 bet1 0 ...0.......................*                 
!  * bet1 al2 bet2 ..0.....................*                  
!  * 0   bet2 al3 bet3 ....................*
!  * .................betn-2 aln-1 betn-1  *
!  * 0 0 0 ..................betn-1 aln    *              
!  * * * * * * * * * * * * * * * * * * * * *
!
!alj = -bj/aj, and betj = sqrt(cj+1/aj/aj+1)
!
!We diagonalize J for Jqj = tqj, with tj being eigenvalues gives grid points
!Now, eigenvectors are used for weights.
! wj = qj[1]**2 * mu0
!
!mu0 = int(a,b) W(x) dx. For GLEG, W(x) = 1, a = -1, b = 1, mu=2.0
! for GLAG, mu = gamma(1+al) for al > -1.
!
!*********************************************************************************************
!*********************************************************************************************
!
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
double precision function aj(j,method)
implicit double precision (a-h,o-z)
integer j
CHARACTER(LEN=*), INTENT(IN) :: method
if (method .eq. 'gleg') then
   aj = (2.d0*j - 1.d0) / j
elseif (method .eq. 'glag') then
   aj = -1.d0 / j
else
   write(*,*) "Provide method either 'gleg' or 'glag'"
endif
return
end


double precision function bj(j,method)
implicit double precision (a-h,o-z)
integer j
CHARACTER(LEN=*), INTENT(IN) :: method
if (method .eq. 'gleg') then
   bj = 0.d0
elseif (method .eq. 'glag') then
   alpha = 1.d0 !change alpha for different value.
   bj = (2.d0*j + alpha - 1.d0) / j
else
   write(*,*) "Provide method either 'gleg' or 'glag'"
endif
return
end


double precision function cj(j,method)
implicit double precision (a-h,o-z)
integer j
CHARACTER(LEN=*), INTENT(IN) :: method
if (method .eq. 'gleg') then
   cj = (j - 1.d0) / j
elseif (method .eq. 'glag') then
   alpha=1.d0 !change alpha for different value
   cj = (j + alpha) / j
else
   write(*,*) "Provide method either 'gleg' or 'glag'"
endif
return
end

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
!Symmetric tridiagonal matrix J.
!
!
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

subroutine GWmat(n,X_mat,method)
implicit none
integer n,i
real(kind=8) X_mat(n,n)
CHARACTER(LEN=*), INTENT(IN) :: method
real(kind=8) aj,bj,cj
external aj,bj,cj
X_mat = 0.d0
if (method .eq. 'gleg') then
   X_mat(1,1) = -1.d0*bj(1,'gleg') / aj(1,'gleg')
   X_mat(1,2) = (cj(2,'gleg') / aj(1,'gleg') / aj(2,'gleg'))**0.5d0
   X_mat(n,n-1) = (cj(n,'gleg') / aj(n-1,'gleg') /   &
&   aj(n,'gleg'))**0.5d0
   X_mat(n,n) = -1.d0*bj(n,'gleg') / aj(n,'gleg')
   do i=2,n-1
    X_mat(i,i) = -1.d0*bj(i,'gleg') / aj(i,'gleg')
    X_mat(i,i+1) = (cj(i+1,'gleg') / aj(i,'gleg') / aj(i+1,'gleg'))**0.5d0
    X_mat(i,i-1) = (cj(i,'gleg') / aj(i-1,'gleg') / aj(i,'gleg'))**0.5d0
   enddo
elseif (method .eq. 'glag') then
   X_mat(1,1) = -1.d0*bj(1,'glag') / aj(1,'glag')
   X_mat(n,n) = -1.d0*bj(n,'glag') / aj(n,'glag')
   X_mat(1,2) = dsqrt(1.d0*cj(2,'glag') / aj(1,'glag') / aj(2,'glag'))
   X_mat(n,n-1) = dsqrt(1.d0*cj(n,'glag') / aj(n-1,'glag') / aj(n,'glag'))
   do i=2,n-1
    X_mat(i,i) = -1.d0*bj(i,'glag') / aj(i,'glag')
    X_mat(i,i+1) = dsqrt(1.d0*cj(i+1,'glag') / aj(i,'glag') / aj(i+1,'glag'))
    X_mat(i,i-1) = dsqrt(1.d0*cj(i,'glag') / aj(i-1,'glag') / aj(i,'glag'))
   enddo
else
    write(*,*) "Provide method either 'gleg' or 'glag'"
endif
return
end

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
!Diagonalization using LAPACK and BLAS. Need lapack and blas library to run this
!code.
!
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

subroutine diagonalize(n,m,xmat,eig,wi,xvec)
implicit double precision (a-h,o-z)
integer n
real(kind=8) eig(n), xvec(n,n),DUMMY(1,1),wi(n),WORK(m),xmat(n,n)
external DGEEV
call DGEEV('N', 'V', n, xmat, n, eig, wi, DUMMY, 1, xvec, n, WORK, m, ok)
if (ok .eq. 0.) then
  write(*,*) "Diagonalization success"
endif
return
end

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
!Calculating grid points and weights for Gauss-quadrature rule.
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine GolWel(n,eig,w,method)
implicit none
integer n,i,m
real(kind=8) w(n),xmat(n,n)
real(kind=8) eig(n),wi(n),xvec(n,n),alpha
CHARACTER(LEN=*), INTENT(IN) :: method

if (method .eq. 'gleg') then
   call GWmat(n,xmat,'gleg')
   m=n*n
   write(*,*) "m: ", m
   call diagonalize(n,m,xmat,eig,wi,xvec)
   do i=1,n
     w(i) = 2.d0*xvec(1,i)**2.d0
   enddo
   
elseif (method .eq. 'glag') then
   call GWmat(n,xmat,'glag')
   call diagonalize(n,m,xmat,eig,wi,xvec)
   alpha=1.d0 !change alpha
   do i=1,n
     w(i) = GAMMA(1+alpha)*xvec(1,i)**2.d0
   enddo
else
   write(*,*) "Provide method either 'gleg' or 'glag'"
endif

return
end

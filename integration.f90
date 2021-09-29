!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!Integration of a function
!Written by Niraj K. Nepal, adapted from Numerical Recipe: Fortran 90
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

program main
implicit double precision (a-h,o-z)
real(kind=8) rslt1, rslt2, rslt3
external func
!call ser(1.d0,2.0d0,1000,f)
!write(*,*) 'array: ',f
call trapezoid(1.d0,2.d0,1000000,rslt1)
call qtrap(1.d0,2.d0,rslt2)
call qgaus(1.d0,2.d0,50,rslt3,'gleg')
write(*,*) 'integration: ', rslt1, rslt2, rslt3
return
end program


! subroutine for creating one-dimensional grid and
! evaluating sinx.
subroutine ser(x1,dx,nstep,fw)
implicit none
real(kind=8) x1,dx
integer nstep,i
real(kind=8), dimension(nstep) :: fw,w
w = 0.d0
fw = 0.d0
do i=1, nstep
  w(i) = w(i) + x1 + dx*(i-1)
  fw(i) = fw(i) + dsin(w(i)) 
end do
return 
end

Double Precision Function func(x)
implicit double precision (a-h,o-z)
func=dsin(x)
return
end


subroutine nrerror(string)
CHARACTER(LEN=*), INTENT(IN) :: String
write(*,*) 'nerror: ',string
STOP 'program terminated by nrerror'
END subroutine nrerror


!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!Trapezoidal rule
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine trapezoid(x1,x2,n,rslt)
implicit double precision (a-h,o-z)
real(kind=8), dimension(n) :: fx
if (n .eq. 1) then
     call ser(x1,x2-x1,1,fx)
     rslt = 0.5d0*(x2-x1)*sum(fx)
else
     dx = (x2 - x1) / n
     call ser(x1,dx,n,fx)
!     write(*,*) size(fx)
     fsum1 = (fx(1) + fx(2))*0.5d0
     fsum2 = sum(fx(2:n-1))
     rslt = dx*(fsum1 + fsum2)
endif
return
end
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! qtrap
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine qtrap(x1,x2,rslt)
implicit none
integer jmax
real(kind=8) x1,x2,rslt,eps,olds
parameter(eps=1d-8,jmax=10000)
integer j
do j=1,jmax
 call trapezoid(x1,x2,j,rslt)
 if (j .gt. 5) then
    if (dabs(rslt-olds) .lt. eps*dabs(olds) .or. &
        & (rslt .lt. 1d-8 .and. olds .lt. 1d-8)) goto 1
 endif
 olds=rslt
enddo
call nrerror('qtrap: too many steps in qtrap')
1 continue
return
end

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! Quadrature rules
!int_(a,b) f(x)dx ~ sum_(i=1,n) wi*f(xi). xi ==> points, wi ==> weight.
!
!int_(a,b) W(x)f(x)dx ~~ sum_(j=1,N) wjf(xj) is exact if f(x) is a polynomial.
! In other words, integrand is smooth enough to approximated by orthogonal
! polynomials Pn(x). Pn(x) has exactly n roots within the interval.
!
!If W(x) = 1, for -1 < x < 1 ==> Gauss-Legendre integration.
!If W(x) = (1 - x**2)**0.5 (-1 < x < 1) ==> Gauss-Chebyshev integration.
!If W(x) = x**alpha (exp(-x)) (0 < x < infty) ==> Gauss-Laguerre
!If W(x) = (1 - x)**alpha (1 + x)**beta (-1 < x < 1) ==> Gauss-Jacobi
!
!The scalar product of 2 functions f and g over a weight function W.
!<f|g> = int_(a,b) W(x)f(x)g(x)dx = 0 if f and g are orthogonal. To find a set
!of polynomials, which is unique for each order j(=0,1,..) (Pj(x)) and are
!orthogonal to each other, we use recurrence relation.
!
!
!
! ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine qgaus(x1,x2,n,rslt,method)
implicit none
real(kind=8) x1,x2,rslt,func
CHARACTER(LEN=*), INTENT(IN) :: method
external func
integer j,n
real(kind=8) dx,xm,xr,w(n),x(n),xt,f
if (method .eq. 'gleg') then
 call gauslegcoeff(-1.d0,1.d0,x,w,n)
endif
xm=0.5d0*(x2+x1)
xr=0.5d0*(x2-x1)
rslt=0.d0
do j=1,n
  dx=xr*x(j)
  xt=xm+dx
  f=func(xt)
  rslt=rslt + w(j)*f
enddo
rslt=xr*rslt
return
end

!*******************************************************************************
!Gauss-Legendre quadrature
!Recurrence relations
!a. (2n-1)xPn-1(x) = nPn(x) + (n-1)Pn-2(x)
!b. (1-x**2)Pn'(x) = nPn-1(x) - nxPn(x)
!c. P-1(x) = 0, P0(x) = 1, n-1==> p2, n-2 ==> p3, n ==> p1
!Points between -1 to +1 and corresponding weight.
!int_(-1,1) f(x)dx ~ sum_(i=1,n) wi*f(xi). xi ==> points, wi ==> weight.
!int_(a,b) f(x)dx ~ sum_(i=1,n) Wi*F(Xi). Xi = (b - a)*xi / 2 + (a + b) / 2
! Wi = (b - a)*wi / 2
!
!For GLEG method, wj = 2 / [(1-x**2)Pn'(x)]. points xi's are the roots of pn(x) 
!
!
!
!
!*******************************************************************************
subroutine gauslegcoeff(x1,x2,x,w,n)
implicit none
integer n
real(kind=8) x1,x2,x(n),w(n)
real(kind=8), parameter :: eps=3.0d-14,pi=3.1415926535897932384626433832795d0
integer i,j,m
real(kind=8) xl,xm,p1,p2,p3,pp,z,z1
m = (n + 1) / 2
xm=0.5d0*(x2+x1)
xl=0.5d0*(x2-x1)
do i=1,m
  z=dcos(pi*(i - 0.25d0) / (n + 0.5d0))
  1 continue
  p1=1.d0
  p2=0.d0
  do j=1,n
   p3=p2
   p2=p1
   p1=((2.d0*j - 1.d0)*z*p2 - (j-1.d0)*p3)/j
  end do
  pp=n*(z*p1 - p2) / (z*z - 1.d0)
  z1=z
  z=z1-p1/pp   !Newton's method
  if (dabs(z-z1) .gt. eps) goto 1
  x(i)=xm-xl*z
  x(n+1-i)=xm+xl*z
  w(i)=2.d0*xl / ((1.d0-z*z)*pp*pp)
  w(n+1-i)=w(i)
enddo
return
end

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!Generating Gauss rule using Golub Welsch method. Any set of polynomials
!{Pj(x)}_(j=1,N) satisfies a 3 term recurrence relation. 
! Pj+1(x) = (aj*x + bj)Pj(x) - cj*Pj-1(x) with P0(x) = 1 and P1(x) = A0*x + B0
! Where aj, bj, and cj are defined as (https://dlmf.nist.gov/18.9)
!
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
! IMPLEMENTATION IN PROGRESS, need diagonalization of matrix.
!*********************************************************************************************
!
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

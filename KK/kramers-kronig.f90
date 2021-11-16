!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!Kramers-kronig relation to obtain real part from the
!imaginary part of the complex function.
! Written by Niraj K. Nepal
!
!For complex function f(x) = Ref(x) + iImf(x)
!
! Ref(x) = (1/pi)*Principal value (PV) of integration_{-inf}^{inf}
! Imf(x')/(x' - x) dx'
!
!
! Imf(x) = -(1/pi)* Principal value (PV) of integration_{-inf}^{inf}
! Ref(x')/(x' - x) dx'
!
! Integration performed using Gauss-Legendre quadrature or Simpson integration.
! For quadrature rule.
! Use integration_qguass.F90 to create grid
! python script to sort the grid. python sorted_grid.py
! cp shorted_grid.dat grid.dat. 
! For accurate results, we need to use specialized logarthmic gaussian
! quadrature. Look at "Numerical evaluation of truncated Kramers-Kronig transformation
! by Frederick W. King.  
!
!
!For linearly spaced frequency interval.
! gfortran -o main kramers-kronig.f90
! python plot.py kk.txt
!
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
program main
implicit double precision (a-h,o-z)
integer, parameter :: nx=10000
real(kind=8) rekk(nx), xt(nx), imf(nx), dfreq, freq(nx), ref(nx)
external aimf
integer ii
imf=0.d0
dfreq=0.001d0
do ii=1,nx
  freq(ii) = (ii-1)*dfreq
  imf(ii) = aimf(freq(ii))
enddo
!call kramerskronig1(aimf,0.d0,10.d0,xt,rekk)
call kramerskronig2(nx,dfreq,freq,imf,ref)
open(1, file='kk.txt', status='new')
do ii=1,nx
! write(1,*) xt(ii), aimf(xt(ii)), rekk(ii), freq(ii), imf(ii), ref(ii)
 write(1,*) freq(ii), imf(ii), ref(ii)
enddo
close(1)
!ccccccccccccccccccccccccccccccccccccccccccccccccccccc
!Reading grid points and weights for the integration.
!ccccccccccccccccccccccccccccccccccccccccccccccccccccc
return
end program

Double Precision Function aimf(y)
implicit double precision (a-h,o-z)
aimf = y*5.d-1 / ((4.d0-y**2.d0)**2.d0 + (5.d-1*y)**2.d0)
return
end

!cccccccccccccccccccccccccccccccccccccccccccccccccccccc
!Calculating real part from the imaginary part
!cccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine kramerskronig1(aimf,x1,x2,yt,rekk)
implicit Double Precision (a-h,o-z)
integer, parameter :: nx = 1000
real(kind=8) x(nx), w(nx), yt(nx),rekk(nx),rslt,a,aimf,y,xm,xr,dx,xt,f
real(kind=8), parameter :: pi=3.1415926535897932384626433832795d0
integer ii,jj,j
external aimf

open(1, file = 'grid.dat', status='old')
do j=1,nx
 read(1,*) x(j), w(j)
enddo
close(1)
rslt=0.d0
rekk=0.d0
xm=0.5d0*(x2+x1)
xr=0.5d0*(x2-x1)
do ii=1,nx
       a = x(ii)*xr + xm
       rslt=0.d0
    do jj=1,nx
       if (jj == ii) CYCLE
       y = x(jj)
       dx=xr*y
       xt=xm+dx
       f = xt*aimf(xt)/(xt**2.d0 - a**2.d0)
       rslt = rslt + w(jj)*xr*f
    enddo
    rekk(ii) = rekk(ii) + (2.d0/pi)*rslt
    yt(ii) = a
enddo
return
end
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!Simpson's integration of order (1/N^4)
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine simpson(nx,dfreq,func,kk)
implicit double precision (a-h,o-z)
integer ii,nx
real(kind=8) dfreq, func(nx), kk
real(kind=8) coef1, coef2, coef3
coef1 = 9.d0/24.d0
coef2 = 28.d0/24.d0
coef3 = 23.d0/24.d0
kk = coef1*func(1) + coef2*func(2) + coef3*func(3)

do ii=4, nx-3
  kk = kk + func(ii)
end do
kk =  kk+coef3*func(nx-2) + coef2*func(nx-1) + coef1*func(nx)
kk = kk * dfreq
return
end
!
!For linear frequency grid.
!
subroutine kramerskronig2(nx,dfreq,freq,imf,ref)
implicit double precision (a-h,o-z)
integer ii, jj
integer nx
real(kind=8), parameter :: pi=3.1415926535897932384626433832795d0
real(kind=8) freq(nx),imf(nx),ref(nx),dfreq,ww,wwp,func(nx),rslt
write(*,*) dfreq
do ii=1,nx
  ww=freq(ii)
  do jj=1,nx
    if (jj == ii) CYCLE
    wwp=freq(jj)
    func(jj) = wwp*imf(jj)/(wwp**2.d0 - ww**2.d0)
  enddo
  call simpson(nx,dfreq,func,rslt)
  ref(ii) = 2.d0/pi * rslt
enddo
return
end

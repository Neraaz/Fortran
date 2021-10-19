!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! Fast Fourier Transform (FFT) to perform Discrete Fourier Transform
! (DFT)
! DFT (N**2) <==> FFT (NLOGN)
! 
! Detail theory: Steve Brunton YouTube lecture.
! https://www.youtube.com/watch?v=E8HeD-MUrjY
!
!
! Download: http://fftw.org/download.html
! Manual: http://manpagez.com/info/fftw3/fftw3-3.3.4/
!
! Installation:
! module load gcc/9.3.0 (If you do not have gcc installed, install it
! first. Simple installation in Unix system.
! ./configure --prefix=/home/tug11655/fortran/fftw3/ --enable-shared
! 
!
! Without --enable-shared, it didn't work.
!
! To enable threaded transformation, add --enable-threads. For OpenMP
! threads, --enable-openmp
!
! For MPI installation, We need to configure using mpi version of
! compiler, link mpi library and header files via LDFLAGS and CPPFLAGS
! variables. We also need to use --enable-mpi.
!
! Look for other configuration options: ./configure --help
!
!Other 2 commands are: (a) make (b) make install
!
!
!Usage:
!module load gcc/9.3.0 (Loading GCC)
!gfortran -L/home/tug11655/fortran/fftw3/lib/ -lfftw3 -I/home/tug11655/fortran/fftw3/include/ fft_program.f90 -o main
!./main
!
!DFT frequency grid: f = [0, 1, ...., n/2 - 1,-n/2, .....-1] / (d*n), if
!n is even and d is sample spacing in real space.
!
!For odd n, [0, 1, ...., (n-1)/2,-(n-1)/2,....-1] / (d*n)
!
!For reference, see numpy.fft.fftfreq() function in python.
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
!Example taken from:
!https://homepage.ntu.edu.tw/~wttsai/fortran/ppt/14.Fast_Fourier_Transform.pdf
!
!
program example_fftw
! Example to call 3-D real FFT routine of FFTW
implicit none
include 'fftw3.f'
integer, parameter :: N=50,L=20,M=30
integer*8 PLAN, PLAN_BAC 
double complex s(L,M,N),e(L,M,N), f(L,M,N),t(L*M*N)
real(kind=8) :: xj, freq(L,M,N),d,freq1(L),freq2(M),freq3(N)
integer j,k,mode,ifftfreq
real(kind=8), parameter :: twopi=2.d0*dacos(-1.d0)
external ifftfreq

write(*,*) "fftw parameters: ", FFTW_FORWARD, FFTW_ESTIMATE, FFTW_BACKWARD
! generating and storing 3D data.
k=L*M*N
do j=1,k
  xj=twopi*real(j)/real(k)
  t(j)=dcos(xj) +0.2d0*dsin(2.d0*xj)
end do
write(*,*) "Original data"
s = reshape(t, (/ L,M,N /))
do j=1,5
 write(*,*) "s: ", s(j,1,1)
enddo
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! Generating and storing 3D data in s(L,M,N) array.
! Forward transform
! Store Fourier coefficients in e(L,M,N) array.
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
call dfftw_plan_dft_3d(PLAN,L,M,N,s,e,FFTW_FORWARD,FFTW_ESTIMATE)
write(*,*) "plan: ", PLAN
call dfftw_execute_dft(PLAN,s,e)
e=e/real(k,KIND=8) ! Normalize
!
!
!!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!!Frequency grid
!!N ==> Number of grid points in real space, d is spacing
!! Now mesh grid of (1 ... L) x (1 ... M), x (1 ... N) in direct space becomes
! FFT mesh grid of (freq1(1) ... freq1(L) x (freq2(1) ... freq2(M)) x (freq3(1)
! ... freq3(N)).
!!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
d=1.d0
do j=1,L
   freq1(j) = ifftfreq(j,L) / (L*d)
enddo
do j=1,M
   freq2(j) = ifftfreq(j,M) / (M*d)
enddo
do j=1,N
   freq3(j) = ifftfreq(j,N) / (N*d)
enddo
 
write(*,*) "Fourier coefficient after forward FFT"
do j=1,5
 write(*,*) "e: ", e(j,1,1)
enddo
!!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!! Backward transform
call dfftw_plan_dft_3d(PLAN_BAC,L,M,N,e,f,FFTW_BACKWARD,FFTW_ESTIMATE) 
write(*,*) "Plan: ", PLAN_BAC
call dfftw_execute_dft(PLAN_BAC,e,f)
write(*,*) "Data after backward FFT"
do j=1,5
  write(*,*) "f: ", f(j,1,1)
end do
!! Destroy the plans
write(*,*) 'ok1'
call dfftw_destroy_plan(PLAN)
write(*,*) 'ok2'
call dfftw_destroy_plan(PLAN_BAC)
write(*,*) 'ok3'
end program example_fftw

integer function ifftfreq(i,n)
implicit none
integer i,n,im
im = i-1
if (im .lt. N/2) then
   ifftfreq = im
else
   ifftfreq = (im - N)
endif
return
end

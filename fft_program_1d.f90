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
! Example to call 1-D real FFT routine of FFTW
implicit none
include 'fftw3.f'
integer, parameter :: N=50
integer*8 PLAN, PLAN_BAC
double complex s(N),e(N), f(N)
real(kind=8) :: xj, freq(N),d
integer j,k,mode,ifftfreq
real(kind=8), parameter :: twopi=2.d0*dacos(-1.d0)
external ifftfreq

write(*,*) "fftw parameters: ", FFTW_FORWARD, FFTW_ESTIMATE, FFTW_BACKWARD
! Discrete data of function f(x)=cos(x)+0.2*sin(2x)
do j=1,N
  xj=twopi*real(j)/real(N)
  s(j)=dcos(xj) +0.2d0*dsin(2.d0*xj)
end do
write(*,*) "Original data"
do j=1,N
  write(*,*) j,s(j)
end do
! Forward transform
call dfftw_plan_dft_1d(PLAN,N,s,e,FFTW_FORWARD,FFTW_ESTIMATE)
write(*,*) "plan: ", PLAN
call dfftw_execute_dft(PLAN,s,e)
e=e/real(N,KIND=8) ! Normalize
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!Frequency grid
!N ==> Number of grid points in real space, d is spacing
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
d=1.d0
do j = 1,N
   freq(j) = ifftfreq(j,N) / (N*d)
enddo
write(*,*) "Fourier coefficient after forward FFT"
do k=1,N
  write(*,*) freq(k),e(k)
end do
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! Backward transform
call dfftw_plan_dft_1d(PLAN_BAC,N,e,f,FFTW_BACKWARD,FFTW_ESTIMATE) 
write(*,*) "Plan: ", PLAN_BAC
call dfftw_execute_dft(PLAN_BAC,e,f)
write(*,*) "Data after backward FFT"
do j=1,N
  write(*,*) j,f(j)
end do
! Destroy the plans
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

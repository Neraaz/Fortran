! using SGEEV to diagonalize Nonsymmetric matrix A.
! taken from https://web.cs.ucdavis.edu/~bai/publications/baidemmeletal06.pdf
!
! For detail about LAPACK functions, please refer to
! http://netlib.org/lapack/explore-html/index.html 
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
program main
implicit none
external SGEEV
real A(4,4),wr(4),wi(4),WORK(16)
real VR(4,4), VL(4,4)
integer i, ok
A = reshape((/ 4., -5., 0., 3., 0., 4., -3., -5., 5., -3., 4., 0., 3., 0., 5., &
& 4. /), (/ 4,4 /) )
!write(*,*), A, shape(A)
call SGEEV('V', 'V', 4, A, 4, wr, wi, VL, 4, VR, 4, WORK, 16, ok)

if (ok .eq. 0) then
  do i=1,4
    write(*,*) wr(i), wi(i), VL(:,i), VR(:,i)
  enddo
else
   write(*,*) "An error occured, ok = ", ok
endif

end program

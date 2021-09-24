!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!solving f(x)=0 for x.
!Written by Niraj K. Nepal adopted from Numerical recipe: Fortran 90
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
program main
implicit double precision (a-h,o-z)
real*8 xb1(100), xb2(100)
logical succes
x1=-1000.d0
x2=1000.d0
!call fbracout(x1,x2,succes,ntry)
!write(*,*) 'x1 and x2', x1,x2,succes,ntry
!if (succes) then
call fbracin(x1,x2,100000000,xb1,xb2,nbb)
write(*,*), nbb
   do i=1,nbb
      sol1=fbisect(xb1(i),xb2(i),1d-12)
      sol2=fflsert(xb1(i),xb2(i),1d-12)
      sol3=fsecant(xb1(i),xb2(i),1d-12)
      sol4=fridder(xb1(i),xb2(i),1d-12)
      write(*,*), 'solution: (bisect,falseposition,secant,ridder)', sol1, sol2, sol3, sol4
   enddo
!else
!    write(*,*) 'range is outside of the solution, either increase ntry in &
!&fbracout or initial guesses'
!endif
return
end program

Double Precision function func(x)
implicit double precision (a-h,o-z)
func = x**2.d0 - 2.d0
return
end

subroutine nrerror(string)
CHARACTER(LEN=*), INTENT(IN) :: String
write(*,*) 'nerror: ',string
STOP 'program terminated by nrerror'
END subroutine nrerror

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!Bracket out and bracket in the roots
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine fbracout(x1,x2,succes,ntry)
implicit double precision (a-h,o-z)
parameter (factor=0.5d0)
external func
logical succes
ntry=0
n=1000
if (x1 .eq. x2) call nrerror('fbracout: you have to guess an initial range')
f1 = func(x1)
f2 = func(x2)
succes=.true.
do j=1, n
  if ((f1 > 0.d0 .and. f2 < 0.d0) .or. (f1 < 0.d0 .and. f2 > 0.d0)) goto 1
  if (dabs(f1) .lt. dabs(f2)) then
     x1=x1+factor*(x1-x2)
     f1=func(x1)
  else
     x2=x2+factor*(x2-x1)
     f2=func(x2)
  end if
  ntry=ntry+1
enddo
succes=.false.
1 continue
return
end


subroutine fbracin(x1,x2,n,xb1,xb2,nbb)
implicit double precision (a-h,o-z)
real*8 xb1(100),xb2(100)
external func
nbb = 0
x = x1
dx = (x2-x1)/n
fp = func(x)
do i=1,n
    x=x+dx
    fc=func(x)
    if(fc*fp .le. 0.d0) then
        nbb=nbb+1
        xb1(nbb)=x-dx
        xb2(nbb)=x
        if (nbb .eq. nb) goto 1
    endif
    fp=fc
!    nb=nbb
enddo
1 continue
return
end

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! root by bisection method
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

Double Precision Function fbisect(x1,x2,xacc)
implicit double precision (a-h,o-z)
external func
maxit = INT(DLOG(dabs(x2-x1)/xacc))
!write(*,*) 'fbisect maxit: ',maxit
fmid=func(x2)
f=func(x1)
if (f*fmid >= 0.d0) call nrerror('fbisect: root is not bracketed')
if (f < 0.d0) then
    fbisect=x1
    dx=x2-x1
else
    fbisect=x2
    dx=x1-x2
endif

do i=1,maxit
    dx=dx*0.5d0
    xmid=fbisect+dx
    fmid=func(xmid)
    if (fmid .le. 0.d0) fbisect=xmid
    if (dabs(dx) < xacc .or. fmid .lt. 1d-12) goto 1
enddo
call nrerror('fbisect: too many bisections')
1 continue
return
end

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!False position root
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
Double Precision Function fflsert(x1,x2,xacc)
implicit double precision (a-h,o-z)
external func
maxit=INT(DLOG(dabs(x2-x1)/xacc))
fl=func(x1)
fh=func(x2)
if ((fl > 0.d0 .and. fh > 0.d0) .or. &
&   (fl < 0.d0 .and. fh < 0.d0)) call &
&   nrerror('fflsert: root must be bracketed between arguments')
if (fl < 0.d0) then
    xl=x1
    xh=x2
else
    xl=x2
    xh=x1
    fl=func(xl)
    fh=func(xh)
endif
dx=xh-xl
do j=1,maxit
    fflsert=xl+dx*fl/(fl-fh)
    f=func(fflsert)
    if (f < 0.d0) then
        del=xl-fflsert
        xl=fflsert
        fl=f
    else
        del=xh-fflsert
        xh=fflsert
        fh=f
    endif
    dx=xh-xl
    if (dabs(del) < xacc .or. f < 1d-12) goto 1
enddo
call nrerror('fflsert: exceed maximum iterations')
1 continue
return
end

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!secant root
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
Double Precision Function fsecant(x1,x2,xacc)
implicit double precision (a-h,o-z)
integer maxit
external func
maxit=100
fl=func(x1)
f=func(x2)
if (dabs(fl) .lt. dabs(f)) then
    fsecant=x1
    xl=x2
    fl=func(xl)
    f=func(fsecant)
else
    xl=x1
    fsecant=x2
endif
do j=1,maxit
   dx=(xl-fsecant)*f/(f-fl)
   xl=fsecant
   fl=f
   fsecant=fsecant+dx
   f=func(fsecant)
   if (dabs(dx) .lt. xacc .or. f < 1d-12) goto 1
enddo
call nrerror('fsecant: exceed maximum iterations')
1 continue
return
end

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!Ridder's method
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
Double Precision Function fridder(x1,x2,xacc)
implicit double precision (a-h,o-z)
external func
maxit=100
fl=func(x1)
fh=func(x2)
if (((fl .gt. 0.d0) .and. (fh .lt. 0.d0)) .or. &
&   (fl .lt. 0.d0 .and. fh .gt. 0.d0)) then
    xl=x1
    xh=x2
    fridder=-1.d30
    do j=1,maxit
      xm=0.5d0*(xl+xh)
      fm=func(xm)
       s=dsqrt(fm**2.d0 - fl*fh)
       if (s < 1d-12) goto 1
       xnew=xm+(xm-xl)*(dsign(1.d0,fl-fh)*fm/s)
       if (dabs(xnew-fridder) .le. xacc) goto 1
       fridder=xnew
       fnew=func(fridder)
       if (fnew .lt. 1d-12) goto 1
       if (dsign(fm,fnew) /= fm) then
         xl=xm
         fl=fm
         xh=fridder
         fh=fnew
       else if (dsign(fl,fnew) /= fl) then
         xh=fridder
         fh=fnew
       else if (dsign(fh,fnew) /= fh) then
         xl=fridder
         fl=fnew
       else
         call nrerror('fridder: never get here')
       endif
       if (dabs(xh-xl) .le. xacc) goto 1
    enddo
    call nrerror('fridder: exceeded max. iterations')
else if (fl .lt. 1d-12) then
    fridder=x1
else if (fh .lt. 1d-12) then
    fridder=x2
else
    call nrerror('fridder: root must be bracketed')
endif
1 continue
return
end


program test_fft
  use fft_mod
  implicit none
  integer, parameter :: n = 64
  integer, parameter :: iter = 1000000
  integer :: ifaxi(100)
  real(8) :: trigxi(100), work(n+1)
  real(8) :: a(n+2)
  integer :: i
  integer(8) :: t1, t2, rate

  do i=1,n
    a(i) = (i-1+14.2D0)**2.D0;
  enddo

  call fftfax_crm(n,ifaxi,trigxi)

  call system_clock(t1)

  do i = 1 , iter
    call fft991_crm(a,work,trigxi,ifaxi,1,n+1,n,1,-1)
    call fft991_crm(a,work,trigxi,ifaxi,1,n+1,n,1,+1)
  enddo

  call system_clock(t2,rate)

  write(*,*) real(t2-t1,8)/real(rate,8)
  write(*,*) a(1)

end program test_fft

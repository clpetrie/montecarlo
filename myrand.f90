program myrand
  use random
  implicit none
  integer, parameter :: num=100000
  integer :: circ
  real :: myrand1, myrand2, pi
  integer :: n

  call init_random_seed()
  circ=0

  do n=1,num
     call random_number(myrand1)
     call random_number(myrand2)
     if (sqrt(myrand1**2+myrand2**2)<1) then
        circ=circ+1
     end if
  end do

  pi=float(circ)/num*4.0
  write(*,*) pi

end program myrand

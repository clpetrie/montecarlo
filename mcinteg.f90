program mcinteg
  use random
  implicit none

  ! set up variables
  integer, parameter:: k2 = selected_real_kind(22)
  real, parameter :: a = 0
  real(k2), parameter :: b = 1E-7
  integer, parameter :: num = 100
  integer, parameter :: steps = 1000
  real, parameter :: h = (b-a)/(num-1)
  real(k2), dimension(num) :: r
  real(k2) :: g, P
  integer :: i
  real(k2) :: pos, myrand, trialpos
  real, dimension(steps) :: steppos
  real :: integral
  real(k2) :: me, hbar, a0, pi, epsilon0, e
  ! declare global variables
  common /variables/ me, hbar, a0, pi, epsilon0, e
  me=9.11E-31
  hbar=1.055E-34
  pi=3.14159
  e=1.602E-19
  epsilon0=8.8542E-12
  a0=4*pi*epsilon0*hbar**2/(me*e**2)
  

  ! make computational grid
  do i=1,num
     r(i)=a+(i-1)*h
  end do

  ! initialize random seed and pick random position
  call init_random_seed
  call random_number(myrand)
  pos=a+(b-a)*myrand

  ! calculate the step positions and do the integral
  integral = 0
  do i=1,steps
     call random_number(myrand)
     trialpos=a+(b-a)*myrand
     if (sin(trialpos)**2 >= sin(pos)**2) then
        pos = trialpos
     else
        call random_number(myrand)
        if (myrand<sin(trialpos)**2/sin(pos)**2) then
           pos=trialpos
        else
           pos=pos
        end if
     end if
     steppos(i)=pos
     integral=integral+g(pos)/P(pos)
  end do

  ! divide final integral by number of sample points
  integral=integral/steps
  integral=integral/e ! convert to eV
  write(*,*) 'E0 = ',integral,' eV'
end program mcinteg

function psi(r)
  implicit none
  integer, parameter:: k2 = selected_real_kind(22)
  real(k2) :: psi, r
  real(k2) :: me, hbar, a0, pi, epsilon0, e
  common /variables/ me, hbar, a0, pi, epsilon0, e
  psi = exp(-r/a0)/sqrt(pi*a0**3)
end function

function psip(r)
  implicit none
  integer, parameter:: k2 = selected_real_kind(22)
  real(k2) :: psip, r
  real(k2) :: me, hbar, a0, pi, epsilon0, e
  common /variables/ me, hbar, a0, pi, epsilon0, e
  psip = -exp(-r/a0)/sqrt(pi*a0**5)
end function

function psipp(r)
  implicit none
  integer, parameter:: k2 = selected_real_kind(22)
  real(k2) :: psipp, r
  real(k2) :: me, hbar, a0, pi, epsilon0, e
  common /variables/ me, hbar, a0, pi, epsilon0, e
  psipp = exp(-r/a0)/sqrt(pi*a0**7)
end function

function g(r)
  implicit none
  integer, parameter:: k2 = selected_real_kind(22)
  real(k2) :: g, r, psi, psip, psipp
  real(k2) :: me, hbar, a0, pi, epsilon0, e
  common /variables/ me, hbar, a0, pi, epsilon0, e
  g=psi(r)*( (-hbar**2/(2*me)) * (psipp(r)+2/r*psip(r)) - e**2/(4*pi*epsilon0*r)*psi(r)) ! this doesn't have the r**2*4*pi because it's only the integrand
end function

function P(r)
  implicit none
  integer, parameter:: k2 = selected_real_kind(22)
  real(k2) :: P, r, psi
  real(k2) :: me, hbar, a0, pi, epsilon0, e
  common /variables/ me, hbar, a0, pi, epsilon0, e
  P=psi(r)**2
end function

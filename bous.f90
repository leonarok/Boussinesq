

module declarations
  implicit none
  private
  double precision,public, allocatable, dimension(:) :: x,x_u,y,y_v
  double precision,public :: xl,yl
  integer,public :: ipref,jpref
  double precision,public, allocatable, dimension(:,:) :: u,v,pc,p,T,rho,mu,gamma,cp
  double precision,public, parameter :: pi=3.1415927
  !double precision,public :: Tamb,radius,perim,h,Tpar
  double precision,public, allocatable, dimension(:,:) :: f_u,f_v
  double precision,public, allocatable, dimension(:,:) :: d_u,d_v
  double precision,public :: m_in,m_out
  real ,public :: relax(6)
  integer,public :: iter,last,npi,npj
end module declarations

module sub
  contains

  subroutine init()
  !initialize variables
  end subroutine

  subroutine grid()
  !define grid
  end subroutine grid

  subroutine fixedbound()
  !specify fixed boundary conditions
  end subroutine fixedbound

  subroutine print()
  !write output
  end subroutine

  subroutine velcorr()
  !correct pressure and velocities
  end subroutine

  subroutine solve()
  !solve linear system with preferred iterative method
  end subroutine

  subroutine ucoeff
  !find coeffiecients for u-equation
  end subroutine

  subroutine vcoeff
  !find coefficients for v-equation
  end subroutine

  subroutine tcoeff
  !find coefficients for T-equation
  end subroutine

  subroutine pcoeff
  !find coeff for pressure correction equation
  end subroutine

end module

program bous
  !main program, solves transient NS
  implicit none
  double precision, allocatable, dimension(:,:)   :: ae,aw,an,as,ap,b
  integer :: k !,it,jt

  call init()
  call grid()
  call fixedbound()

  allocate(ae(npi,npj),aw(npi,npj),an(npi,npj),as(npi,npj),ap(npi,npj),b(npi,npj))
  allocate(f_u(npi,npj),f_v(npi,npj))

  iterlast=1000

  do n = 1,tend/dt
    !calculate coeffiecients
    do iter=1,iterlast
      !solve the linear systems, pressure correction etc.
    end do
  end do

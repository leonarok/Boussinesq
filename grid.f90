subroutine grid()
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  use declarations
  implicit none
  integer :: i,j
  double precision :: dx,dy
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!**** purpose: defining the grid 
!**** see fig. 6.2-6.4 in ref. 1 
!
!
  allocate(x(npi),x_u(npi))
  allocate(y(npj),y_v(npj))
!
!**** purpose: definding the geometrical variables 
!**** see fig. 6.2-6.4 in ref. 1. 
!
!---- length of the area in the x- and y direction 
!---- in milimeter
  xl=500
  yl=500
!
!---- length of volume element
!
  dx=xl/real(npi-2)
  dy=yl/real(npj-2)
!
!
!---- length variable for the scalar points in the x direction
!
  x(1)=0.
  x(2)=0.5*dx
  do i=3,npi-1
    x(i)=x(i-1)+dx
  end do
  x(npi)=x(npi-1)+0.5*dx
!
!---- length variable for the scalar points fi(i,j) in the y direction
!
  y(1)=0.
  y(2)=0.5*dy
  do j=3,npj-1
    y(j)=y(j-1)+dy
  end do
  y(npj)=y(npj-1)+0.5*dy
!
!---- length variable for the velocity components u(i,j) in the x direction
!
  x_u(1)=0.
  x_u(2)=0.
  do i=3,npi
    x_u(i)=x_u(i-1)+dx
  end do
!
!---- length variable for the velocity components v(i,j) in the y direction
!
  y_v(1)=0.
  y_v(2)=0.
  do j=3,npj
    y_v(j)=y_v(j-1)+dy
  end do   
!
end subroutine grid
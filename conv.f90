subroutine conv()
  use declarations
  implicit none
  integer :: i,j

  do i=2, npi
  	do j=2, npi
      f_u(i,j)=(rho*(x(i)-x_u(i))+ &
        rho∗(x_u(i)-x(i-1)))∗u(i,j)/(x(i)-x(i-1))
      f_v(i,j)=(rho∗(y(j)-y_v(j))+ &
        rho∗(y_v(j)-y(j-1)))∗v(i,j)/(y(j)-y(j-1))!=f(i,j)
    end do
  end do
end subroutine conv

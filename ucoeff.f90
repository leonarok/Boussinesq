subroutine ucoeff(ae, aw, an, as, ap, b)
use declarations 
  implicit none 
  double precision, dimension(:,:), intent(out) :: ae, aw, an, as, ap, b
  integer :: i,j
  real :: fw, fe, fs, fn, dw, de, ds, dn, areaw, areae, areas, arean, sp, su

  call conv()

    do i=, npi-1 
     do j=2, npj-1
       !

     areaw = y_v ( j +1) -y_v ( j ) ! s e e f i g . 6 . 3
     areae = areaw
     areas = x(i)-x(i-1)
     arean = areas      

     fw=((f_u(i,j)+f_u(i-1,j))/2)∗areaw
     fe=((f_u(i+1,j)+f_u(i,j))/2)∗areae
     fs=((f_v(i,j)+f_v(i-1,j))/2)∗areas
     fn=((f_v(i,j+1)+f_v(i-1,j+1))/2)∗arean      

         

     dw=(mu/(x_u(i)-x_u(i-1)))∗areaw
     de=(mu/(x_u(i+1)-x_u(i)))∗areae
     ds=(mu/((y(j)-y(j-1))))∗areas
     dn=(mu/((y(j+1)-y(j))))∗arean      

     sp=0.
     su=0.      

     aw(i,j)=dw+max(fw,0.)
     ae(i,j)=de+max(0.,-fe)
     as(i,j)=ds+max(fs,0.)
     an(i,j)=dn+max(0.,-fn)      

     ap(i,j)=aw(i,j)+ae(i,j)+as(i,j)+an(i,j)+fe-fw+fn-fs-sp      

     d_u(i,j)=areaw∗relax(1)/ap(i,j)      

     b(i,j)=(p(i-1,j)-p(i,j))∗areaw+su      

     ap(i,j)=ap(i,j)/relax(1)
     
     b(i,j)=b(i,j)+(1-relax(1))∗ap(i,j)∗u(i,j)      

    end do
  end do
end subroutine ucoeff
subroutine pccoeff(ae,aw,an,as,ap,b)

use declarations
implicit none
double precision, dimension(:,:),intent(out)::ae,aw,an,as,ap,b
integer::i,j
real :: areaw,areae,areas,arean,sp

call conv()

do i=2,npi-1
  do j=2,npj-1

    areaw=y_v(j+1)-y_v(j)
    areae=areaw
    areas=x_u(i+1)-x_u(i)
    arean=areas   

    b(i,j)=f_u(i,j)∗areaw-f_u(i+1,j)∗areae+   &
      f_v(i,j)∗areas-f_v(i,j+1)∗arean
    sp=0.   

    ae(i,j)=rho*d_u(i+1,j)*areae
    aw(i,j)=rho*d_u(i,j)*areaw
    an(i,j)=rho*d_v(i+1,j)*arean
    as(i,j)=rho*d_v(i+1,j)*areas

    ap(i,j)=ae(i,j)+aw(i,j)+an(i,j)+as(i,j)-sp

    pc(i,j)=0

  enddo
enddo

end subroutine pccoeff
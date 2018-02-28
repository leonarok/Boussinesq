subroutine tcoeff(ae,aw,an,as,ap,b)

use declarations
implicit none
double precision,dimension(:,:),intent(out)::ae,aw,an,as,ap,b
integer::i,j
real::fw,fe,fs,fn,dw,de,ds,dn,areaw,areae,areas,arean,sp,su

call conv()

  do i=2,npi-1
    do j=2,npj-1

      areaw=y_v(j+1)-y_v(j)
      areae=areaw
      areas=x_u(i+1)-x_u(i)
      arean=areas     

      fw=f_u(i,j)∗cp∗areaw
      fe=f_u(i+1,j)∗cp∗areae
      fs=f_v(i,j)∗cp∗areas
      fn=f_v(i,j+1)∗cp∗arean      

      dw=((gamma)/(gamma∗   &
        (x(i)-x_u(i))+gamma∗(x_u(i)-x(i-1))))∗areaw
      de=((gamma∗gamma)/(gamma∗&
        (x(i+1)-x_u(i+1))+gamma∗(x_u(i+1)-x(i))))∗areae
      ds=((gamma∗gamma)/(gamma∗ &
        (y(j)-y_v(j))+gamma∗(y_v(j)-y(j-1))))∗areas
      dn=((gamma∗gamma)/(gamma∗   &
        (y(j+1)-y_v(j+1))+gamma∗(y_v(j+1)-y(j))))∗arean     

      sp=0.
      su=0.     

      aw(i,j)=dw+max(fw,0.)
      ae(i,j)=de+max(0.,-fe)
      as(i,j)=ds+max(fs,0.)
      an(i,j)=dn+max(0.,-fn)      

      ap(i,j)=aw(i,j)+ae(i,j)+as(i,j)+an(i,j)+fe-fw+fn-fs-sp      

      b(i,j)=su     

      enddo
    enddo     

end subroutine tcoeff
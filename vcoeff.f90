subroutine vcoeff(ae,aw,an,as,ap,b)
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
use declarations
implicit none
doubleprecision,dimension(:,:),intent(out)::ae,aw,an,as,ap,b
integer::i,j
real::fw,fe,fs,fn,dw,de,ds,dn,areaw,areae,areas,arean,sp,su

call conv()
!
do i=2,npi-1
	do j=3,npj-1
    areaw=y(j)-y(j-1)
    areae=areaw
    areas=x_u(i+1)-x_u(i)
    arean=areas    
    
    fw=((f_u(i,j)+f_u(i,j-1))/2)∗areaw
    fe=((f_u(i+1,j)+f_u(i+1,j-1))/2)∗areae
    fs=((f_v(i,j)+f_v(i,j-1))/2)∗areas
    fn=((f_v(i,j)+f_v(i,j+1))/2)∗arean    

    dw=(mu/((x(i)-x(i-1))))∗areaw
    de=(mu/((x(i+1)-x(i))))∗areae
    ds=(mu/(y_v(j)-y_v(j-1)))∗areas
    dn=(mu/(y_v(j+1)-y_v(j)))∗arean    

    sp=0.
    su=0.    

    aw(i,j)=dw+max(fw,0.)
    ae(i,j)=de+max(0.,-fe)
    as(i,j)=ds+max(fs,0.)
    an(i,j)=dn+max(0.,-fn)    

    ap(i,j)=aw(i,j)+ae(i,j)+as(i,j)+an(i,j)+fe-fw+fn-fs-sp    

    d_v(i,j)=areas∗relax(2)/ap(i,j)    

    b(i,j)=(p(i,j-1)-p(i,j))∗areas+su    

    ap(i,j)=ap(i,j)/relax(2)

    b(i,j)=b(i,j)+(1-relax(2))∗ap(i,j)∗v(i,j)    

    enddo
	enddo

end subroutine vcoeff
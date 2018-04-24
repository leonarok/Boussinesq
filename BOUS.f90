!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!-----------------------------------------------------------------------------!
!-------------Module containing central variables for the program-------------!
!-----------------------------------------------------------------------------!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module declarations

	implicit none
	private
	double precision,public, allocatable, dimension(:) :: x,x_u,y,y_v,er
	double precision,public :: xl,yl
	integer,public :: ipref,jpref
	double precision,public,allocatable,dimension(:,:) :: u,v,pc,p,T
	double precision,public,allocatable,dimension(:,:) :: u_tent,v_tent,T_tent
	double precision,public, parameter :: pi=3.1415927
	double precision,public :: radius,perim,h,Tpar,rho,mu,nu,nu_a,k,k_a,cp
	double precision,public, allocatable, dimension(:,:) :: f_u,f_v
	double precision,public, allocatable, dimension(:,:) :: d_u,d_v
	double precision,public :: dx,dy
	real ,public :: relax(6)
	integer,public :: iter,last,npi,npj,nsteps
	double precision,public :: cup_width,k_cup,T_amb,g_const,alpha_l,tend,dt
	double precision,public :: h_coVERT,h_coHORI,h_ci,beta_l,beta_a,h_co
	double precision,public :: Ra,Pr_l,Pr_a
end module declarations

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!-----------------------------------------------------------------------------!
!----------------Module containing subroutines for the program----------------!
!-----------------------------------------------------------------------------!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module sub

	contains
	
	SUBROUTINE tick(t)
        INTEGER, INTENT(OUT) :: t

        CALL system_clock(t)
    END SUBROUTINE tick

    ! returns time in seconds from now to time described by t
    REAL FUNCTION tock(t)
        INTEGER, INTENT(IN) :: t
        INTEGER :: now, clock_rate

        call system_clock(now,clock_rate)

        tock = real(now - t)/real(clock_rate)
    END FUNCTION tock

	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!-------------------------------------------------------------------------!
	!--------Initializes all parameters and sets number of grid points--------!
	!-------------------------------------------------------------------------!
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	subroutine init()
	
		!----------------!
		! Initialization !
		!----------------!
		use declarations
		implicit none
		
		!-----------------------------------------------------------------!
		! Set number of grid points in each direction and allocate arrays !
		!-----------------------------------------------------------------!
		npi=92
	  	npj=122
	  	last=400	
		allocate(u(npi,npj),v(npi,npj),p(npi,npj),pc(npi,npj))
		allocate(T(npi,npj))
		allocate(u_tent(npi,npj),v_tent(npi,npj),T_tent(npi,npj))
		allocate(d_u(npi,npj),d_v(npi,npj),er(last))
	  
		!--------------------------------------!
		! Set timestep and simulation end time !
		!--------------------------------------!
		tend=3600
		dt=10
		nsteps=ceiling(tend/dt)

		!---------------------------------!
		! Set reference for zero pressure !
		!---------------------------------!
		ipref=npi/2
		jpref=2
	
		!----------------------------------!
		! Initialize simulation parameters !
		!----------------------------------!
		T_amb=273.16+20
		u=0.00
		u_tent=0.00
		v=0.00
		v_tent=0.00
		p=0.00 !initial pressure estimate
		pc=0.00
		T(1:npi,1:npj)=273.16+80.00
		d_u=0.
		d_v=0.
	
		!---------------------!
		! Set constant values !
		!---------------------!
	 	g_const=9.81			! Standard gravity
	 	alpha_l=0.000069			! Thermal expansion coeff. of liquid
	 	rho=1000	! Density
	 	mu=0.001   ! Dynamic viscosity
	 	nu=mu/rho	! Kinematic viscosity
	 	k=0.64  	! Thermal conductivety 
	 	cp=4181   	! Specific heat of liquid
	 	beta_l=k/(rho*cp)	 	! Thermal diffusivity of liquid
	 	Pr_l=cp*mu/k 			! Prandtl number for liquid

	  
	  	!-----------------------------------!
		! Set cup/wall heat loss parameters !
		!-----------------------------------!
	    cup_width=0.008
	    k_cup=2.09
	   	Pr_a=0.75				! Prandtl number for air
	  	beta_a=0.000019			! Thermal diffusivity of air
	  	k_a=0.025
	  	nu_a=0.0000148
	  	
	  	!---------------------------!
		! Set relaxation parameters !
		!---------------------------!
		relax(1)=0.85 	! Relax u
		relax(2)=0.85 	! Relax v
		relax(4)=0.25 	! Relax p
		relax(5)=1.0 	! Relax T
		
	end subroutine init
	
	
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!-------------------------------------------------------------------------!
	!-----------------Defines grid and computational domain-------------------!
	!-------------------------------------------------------------------------!
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	subroutine grid()
		
		!----------------!
		! Initialization !
		!----------------!
		use declarations
	  	implicit none
	 	integer :: i,j
		
		!-------------------------!
		! Allocate length vectors !
		!-------------------------!
	  	allocate(x(npi),x_u(npi))
	  	allocate(y(npj),y_v(npj))

		!----------------------------------!
		! Set size of computational domain !
		!----------------------------------!
	  	xl=0.08-cup_width*2
	  	yl=0.09-cup_width
		
		!---------------------------!
		! Calculate inner grid size !
		!---------------------------!
	  	dx=xl/real(npi-2)
	  	dy=yl/real(npj-2)

		!------------------------------------------------------!
		! Length variable for scalar points in the x-direction !
		!------------------------------------------------------!
	  	x(1)=0.
	  	x(2)=0.5*dx
	  	do i=3,npi-1
	    	x(i)=x(i-1)+dx
	  	end do
	  	x(npi)=x(npi-1)+0.5*dx

		!------------------------------------------------------!
		! Length variable for scalar points in the y-direction !
		!------------------------------------------------------!
	  	y(1)=0.
	  	y(2)=0.5*dy
	  	do j=3,npj-1
	    	y(j)=y(j-1)+dy
	  	end do
	  	y(npj)=y(npj-1)+0.5*dy

		!-------------------------------------------------------------------!
		! Length variable for velocity components u(i,j) in the x-direction !
		!-------------------------------------------------------------------!
	  	x_u(1)=0.
	  	x_u(2)=0.
	  	do i=3,npi
	    	x_u(i)=x_u(i-1)+dx
	  	end do
		
		!-------------------------------------------------------------------!
		! Length variable for velocity components v(i,j) in the x-direction !
		!-------------------------------------------------------------------!
	 	y_v(1)=0.
	  	y_v(2)=0.
	  	do j=3,npj
	    	y_v(j)=y_v(j-1)+dy
	  	end do   
	  	
	end subroutine grid
	      
	
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!-------------------------------------------------------------------------!
	!--------Specifies fixed boundary values that remains untouched-----------!
	!-------------------------------------------------------------------------!
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	subroutine fixedbound()
	
		!----------------!
		! Initialization !
		!----------------!
		use declarations
		implicit none

		!--------------------------------------!
		! Impermeability on east and west wall !
		!--------------------------------------!
	  	u(1,:)=0.
	  	u(2,:)=0.
	  	u(npi,:)=0.
	  	
	  	!----------------------------------------!
		! Impermeability on north and south wall !
		!----------------------------------------!
	  	v(:,1)=0.
	  	v(:,2)=0.
	  	v(:,npj)=0.

		!---------------------------------!
		! No-slip on north and south wall !
		!---------------------------------!
	  	u(:,1)=0.
	  	u(:,npj)=0.
	
		!-------------------------------!
		! No-slip on east and west wall !
		!-------------------------------!
	 	v(1,:)=0.
	  	v(npi,:)=0.
	  
	end subroutine fixedbound
	
	
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!-------------------------------------------------------------------------!
	!-------------Specifies boundary values for every iteration---------------!
	!-------------------------------------------------------------------------!
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	subroutine bound()
	 	
	 	!----------------!
		! Initialization !
		!----------------!
	    use declarations
	    implicit none
	    		
	  	!------------------------------------------------------!
		! von Neumann BC for temperature on east and west wall !
		!------------------------------------------------------!
	    T_tent(npi,2:npj-1)=T_tent(npi-1,2:npj-1)
	    T_tent(1,2:npj-1)=T_tent(2,2:npj-1)
	
	  	!--------------------------------------------------------!
		! von Neumann BC for temperature on north and south wall !
		!--------------------------------------------------------!
	  	T_tent(2:npi-1,1)=T_tent(2:npi-1,2)
	  	T_tent(2:npi-1,npj)=T_tent(2:npi-1,npj-1)
	
	end subroutine bound
	
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!-------------------------------------------------------------------------!
	!----------Corrects pressure and velocities by SIMPLE algotihm------------!
	!-------------------------------------------------------------------------!
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	subroutine velcorr()
		
		!----------------!
		! Initialization !
		!----------------!
		use declarations
		implicit none
	  	integer :: i,j

	 	do i=2,npi-1
	  		do j=2,npj-1
	  		
				!-------------------------------------!
				! Pressure correction with relaxation !
				!-------------------------------------!
	      		p(i,j)=p(i,j)+relax(4)*(pc(i,j)-pc(ipref,jpref))
				
				!---------------------!
				! Velocity correction !
				!---------------------!
	      		if(i.ne.2) u_tent(i,j)=u_tent(i,j)+d_u(i,j)*(pc(i-1,j)-pc(i,j))
	      		if(j.ne.2) v_tent(i,j)=v_tent(i,j)+d_v(i,j)*(pc(i,j-1)-pc(i,j))
	      
	    	end do
	  	end do
	  	
	  	!------------------------------------!
		! Extrapolate pressure at boundaries !
		!------------------------------------!   
	  	p(1,:) = p(2,:)  
	  	p(:,1) = p(:,2) 
	  	p(npi,:) = p(npi-1,:)
	  	p(:,npj) = p(:,npj-1)
	  	     
	end subroutine velcorr
	
	
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!-------------------------------------------------------------------------!
	!-----Prints out results to output directory for the current timestep-----!
	!-------------------------------------------------------------------------!
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	subroutine print(time)
	
		!----------------!
		! Initialization !
		!----------------!
		use declarations
		implicit none
		integer :: i,j
		real :: time
		character(len=30) :: rowfmt
		character(len=30) :: filename
		
		!------------------------!
		! Establish write format !
		!------------------------!
		write(rowfmt,'(A,I4,A)') '(',npj,'(2X,E15.6E2))'
		
		!-------------------!
		! Open output files !
		!-------------------!
		write(filename,'("output/u/u_",f8.2,".dat")') time
	    open(100,file=filename,status="unknown",RECL=(17*npj+120))
	    
	    write(filename,'("output/v/v_",f8.2,".dat")') time
	    open(101,file=filename,status='unknown',RECL=(17*npi+120))
	    
	    write(filename,'("output/temp/temp_",f8.2,".dat")') time
	   	open(12,file=filename,status='unknown',RECL=(17*npj+120))
	   	
	   	write(filename,'("output/p/p_",f8.2,".dat")') time
	    open(13,file=filename,status='unknown',RECL=(17*npj+120))

	    !----------------------------------------!
		! Write x- and y-vectors to output files !
		!----------------------------------------!
	    open(110,file='output/x.dat',status='unknown')    
	    open(111,file='output/y.dat',status='unknown')
	    write(110,'(1600F14.7)') x(:)
	    write(111,'(1600F14.7)') y(:)
	
	
		!---------------------------------!
		! Write u-velocity to output file !
		!---------------------------------!
	    write(100,FMT=rowfmt) u(2,:)
	    do i=2,npi-1
	      	write(100,FMT=rowfmt) (u(i,:)+u(i+1,:))/2.
	    end do
	    write(100,FMT=rowfmt) u(npi,:)
		
		!---------------------------------!
		! Write v-velocity to output file !
		!---------------------------------!
	    write(101,FMT=rowfmt) v(:,2)
	    do j=2,npj-1
	      	write(101,FMT=rowfmt) (v(:,j)+v(:,j+1))/2.
	    end do
	    write(101,FMT=rowfmt) v(:,npj)

		!----------------------------------!
		! Write temperature to output file !
		!----------------------------------!
	    do i=1,npi
	    	write(12,FMT=rowfmt) T(i,:)
	    end do
		
		!-------------------------------!
		! Write pressure to output file !
		!-------------------------------!
	    do i=1,npi
	      	write(13,FMT=rowfmt) p(i,:)
	    end do

		!-------------!
		! Close files !
		!-------------!
	    close(100)
	    close(101)
	    close(110)
	    close(111)
	    close(12)
	    close(13)
	    
	end subroutine print
	
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!-------------------------------------------------------------------------!
	!--------Solves linear system of eq.s by line Gauss-Seidel method---------!
	!-------------------------------------------------------------------------!
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	subroutine solve(fi,b,ae,aw,an,as,ap,istart,iend,jstart,jend)
		
		!----------------!
		! Initialization !
		!----------------!
	  	implicit none
	  	double precision, dimension(:,:), intent(in out) :: fi
	  	double precision, dimension(:,:), intent(in) :: b
	  	double precision, dimension(:,:), intent(in) :: ae,aw,an,as,ap
	  	integer, intent(in) :: istart,iend,jstart,jend
	  	integer :: i,j
	  	double precision, allocatable, dimension(:) :: ath, cmth
	  	double precision :: cth
		
		!-----------------------!
		! Allocate line vectors !
		!-----------------------!	
	  	allocate(ath(iend),cmth(iend))
	  	
		!----------------------------------------!
		! Solving the (e-w) lines from the south !
		!----------------------------------------!
	  	do j=jstart+1,jend-1 
			!-----------------------------!
			! Set values at west boundary !
			!-----------------------------!
	    	ath(istart)=0.   
	    	cmth(istart)=fi(istart,j) 
	    	
	    	!----------------------!
			! Forward substitution !
			!----------------------!
	    	do i=istart+1,iend-1
	      		ath(i)=ae(i,j)/(ap(i,j)-aw(i,j)*ath(i-1))
	      		cth=an(i,j)*fi(i,j+1)+as(i,j)*fi(i,j-1)+b(i,j)   
	      		cmth(i)=(aw(i,j)*cmth(i-1)+cth)/(ap(i,j)-aw(i,j)*ath(i-1))
	    	end do   
	    	
	    	!-----------------------!
			! Backward substitution !
			!-----------------------!
	    	do i=iend-1,istart+1,-1
	      		fi(i,j)=ath(i)*fi(i+1,j)+cmth(i)
	    	end do
	  	end do

		!----------------------------------------!
		! Solving the (e-w) lines from the north !
		!----------------------------------------!
	  	do j=jend-2,jstart+1,-1 
			
			!-----------------------------!
			! Set values at west boundary !
			!-----------------------------!
	    	ath(istart)=0. 
	    	cmth(istart)=fi(istart,j)
			
			!----------------------!
			! Forward substitution !
			!----------------------!
	    	do i=istart+1,iend-1
	      		ath(i)=ae(i,j)/(ap(i,j)-aw(i,j)*ath(i-1))
	      		cth=an(i,j)*fi(i,j+1)+as(i,j)*fi(i,j-1)+b(i,j)  
	      		cmth(i)=(aw(i,j)*cmth(i-1)+cth)/(ap(i,j)-aw(i,j)*ath(i-1))
	    	end do   

			!-----------------------!
			! Backward substitution !
			!-----------------------!
	    	do i=iend-1,istart+1,-1  
	      		fi(i,j)=ath(i)*fi(i+1,j)+cmth(i)
	    	end do
	  	end do      
		
		!-------------------------!
		! Deallocate line vectors !
		!-------------------------!
	  	deallocate(ath,cmth)
		
		!-----------------------!
		! Allocate line vectors !
		!-----------------------!
	  	allocate(ath(jend),cmth(jend))
		
		!---------------------------------------!
		! Solving the (n-s) lines from the west !
		!---------------------------------------!
	  	do i=istart+1,iend-1
			
			!------------------------------!
			! Set values at south boundary !
			!------------------------------!
	    	ath(jstart)=0.   
	    	cmth(jstart)=fi(i,jstart) 
	    	
	    	!----------------------!
			! Forward substitution !
			!----------------------!
	    	do j=jstart+1,jend-1
	      		ath(j)=an(i,j)/(ap(i,j)-as(i,j)*ath(j-1))
	      		cth=ae(i,j)*fi(i+1,j)+aw(i,j)*fi(i-1,j)+b(i,j)  
	      		cmth(j)=(as(i,j)*cmth(j-1)+cth)/(ap(i,j)-as(i,j)*ath(j-1))
	    	end do   
			
			!-----------------------!
			! Backward substitution !
			!-----------------------!
	    	do j=jend-1,jstart+1,-1
	      		fi(i,j)=ath(j)*fi(i,j+1)+cmth(j)
	    	end do   
	  	end do

		!---------------------------------------!
		! Solving the (n-s) lines from the east !
		!---------------------------------------!
	  	do i=iend-2,istart+1,-1
			
			!------------------------------!
			! Set values at south boundary !
			!------------------------------!
	    	ath(jstart)=0. 
	    	cmth(jstart)=fi(i,jstart)
   			
   			!----------------------!
			! Forward substitution !
			!----------------------!      
	    	do j=jstart+1,jend-1
	      		ath(j)=an(i,j)/(ap(i,j)-as(i,j)*ath(j-1))
	      		cth=ae(i,j)*fi(i+1,j)+aw(i,j)*fi(i-1,j)+b(i,j)  
	      		cmth(j)=(as(i,j)*cmth(j-1)+cth)/(ap(i,j)-as(i,j)*ath(j-1))
	    	end do   

			!-----------------------!
			! Backward substitution !
			!-----------------------!
	    	do j=jend-1,jstart+1,-1
	      		fi(i,j)=ath(j)*fi(i,j+1)+cmth(j) 
	    	end do   
	  	end do  
	  	
	  	!-------------------------!
		! Deallocate line vectors !
		!-------------------------!
	  	deallocate(ath,cmth)
	  	
	end subroutine solve
	
	
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!-------------------------------------------------------------------------!
	!--------------Calculates coefficients for the u-equation-----------------!
	!-------------------------------------------------------------------------!
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	subroutine ucoeff(ae,aw,an,as,ap,b)

		!----------------!
		! Initialization !
		!----------------!
		use declarations
	  	implicit none
	  	double precision, dimension(:,:), intent(out) :: ae,aw,an,as,ap,b
	  	integer :: i,j
	  	real :: fw,fe,fs,fn,dw,de,ds,dn,areaw,areae,areas,arean
	  	real :: sp,su,ap_zero,ap_tilde
	  	real :: r_w,r_e,dx_cell,dy_cell,volume
		
		!-----------------------------!
		! Calculates convective terms !
		!-----------------------------!
	  	call conv()
	  
	  	!-----------------------!
		! Iterate through array !
		!-----------------------!
	  	do i=3,npi-1
	    	do j=2,npj-1

				!-------------------------!
				! Find area of cell faces !
				!-------------------------!
!	      		areaw=y_v(j+1)-y_v(j) 
!	      		areae=areaw
!	      		areas=x(i)-x(i-1)
!	      		arean=areas
				dx_cell=y_v(j+1)-y_v(j)
				dy_cell=x(i)-x(i-1)
				r_w=abs(xl/2-x(i-1))
				r_e=abs(xl/2-x(i))
				
				areaw=pi*r_w*dy_cell
	      		areae=pi*r_e*dy_cell
	      		areas=pi/2*abs(r_w**2-r_e**2)
	      		arean=areas
	      		volume=areas*dy
			

				!---------------------------!
				! Find convective mass flux !
				!---------------------------!
	        	fw=0.5*( f_u(i,j)+f_u(i-1,j) ) * areaw
	        	fe=0.5*( f_u(i+1,j)+f_u(i,j) ) * areae
	        	fs=0.5*( f_v(i,j)+f_v(i-1,j) ) * areas
	        	fn=0.5*( f_v(i,j+1)+f_v(i-1,j+1) ) * arean

	      		!----------------------------!
				! Find diffusion conductance !
				!----------------------------!
	        	dw=(mu/(x_u(i)-x_u(i-1)))*areaw
	        	de=(mu/(x_u(i+1)-x_u(i)))*areae
	        	ds=(mu/(y(j)-y(j-1)))*areas
	        	dn=(mu/(y(j+1)-y(j)))*arean
	
	      		!-------------------!
				! Find source terms !
				!-------------------!
	        	sp=0.
	        	su=0.

	      		!----------------------------------!
				! Establish neighbour coefficients !
				!----------------------------------!
	        	aw(i,j)=dw+fw/2
	        	ae(i,j)=de-fe/2
	        	as(i,j)=ds+fs/2
	        	an(i,j)=dn-fn/2
	    
	      		!------------------------------!
				! Establish center coefficient !
				!------------------------------!
	        	ap_zero= rho*volume/dt
	        	ap_tilde=aw(i,j)+ae(i,j)+as(i,j)+an(i,j)+fe-fw+fn-fs-sp
	        	ap(i,j)=ap_zero+ap_tilde
	
	      		!-----------------------------------------------!
				! Calculates coefficient term for use in pc eq. !
				!-----------------------------------------------!
	        	d_u(i,j)=(areaw+areae)/2*relax(1)/ap(i,j)
	
				!-----------------------------!
				! Establish total source term !
				!-----------------------------!
	        	b(i,j)=(p(i-1,j)-p(i,j))*(areaw+areae)/2+ap_zero*u(i,j)+su
	
	      		!--------------------------------------------------------------!
				! Introducing relaxation and adjusting source term accordingly !
				!--------------------------------------------------------------!
	        	ap(i,j)=ap(i,j)/relax(1)
	        	b(i,j)=b(i,j)+(1-relax(1))*ap(i,j)*u_tent(i,j)
	    end do
	  end do   
	  
	end subroutine ucoeff
	
	
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!-------------------------------------------------------------------------!
	!--------------Calculates coefficients for the v-equation-----------------!
	!-------------------------------------------------------------------------!
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	subroutine vcoeff(ae,aw,an,as,ap,b)
		
		!----------------!
		! Initialization !
		!----------------!
	   	use declarations
	    implicit none
	    double precision, dimension(:,:), intent(out) :: ae,aw,an,as,ap,b
	    integer :: i,j
	    real :: fw,fe,fs,fn,dw,de,ds,dn,areaw,areae,areas,arean,sp,su
	    real :: ap_tilde,ap_zero
	    real :: r_w,r_e,dx_cell,dy_cell,volume
	    
	    !-----------------------------!
		! Calculates convective terms !
		!-----------------------------!
	    call conv()
	
		!-----------------------!
		! Iterate through array !
		!-----------------------!
	    do i=2,npi-1
	      	do j=3,npj-1
	      	
	      		!-------------------------!
				! Find area of cell faces !
				!-------------------------!
!	        	areaw=y(j)-y(j-1)
!	        	areae=areaw
!	        	areas=x_u(i+1)-x_u(i)
!	        	arean=areas
				dx_cell=y(j)-y(j-1)
				dy_cell=x_u(i+1)-x_u(i)
				r_w=abs(xl/2-x_u(i))
				r_e=abs(xl/2-x_u(i+1))
				
				areaw=pi*r_w*dy_cell
	      		areae=pi*r_e*dy_cell
	      		areas=pi/2*abs(r_w**2-r_e**2)
	      		arean=areas
	      		volume=areas*dy
	
	      		!---------------------------!
				! Find convective mass flux !
				!---------------------------!
	        	fw=0.5*( f_u(i,j)+f_u(i,j-1) ) * areaw
	        	fe=0.5*( f_u(i+1,j)+f_u(i+1,j-1) ) * areae
	        	fs=0.5*( f_v(i,j)+f_v(i,j-1) ) * areas
	        	fn=0.5*( f_v(i,j)+f_v(i,j+1) ) * arean
	
	     		!----------------------------!
				! Find diffusion conductance !
				!----------------------------!
	        	dw=(mu/(x(i)-x(i-1)))*areaw
	        	de=(mu/(x(i+1)-x(i)))*areae
	        	ds=(mu/(y_v(j)-y_v(j-1)))*areas
	        	dn=(mu/(y_v(j+1)-y_v(j)))*arean
	
	      		!-------------------!
				! Find source terms !
				!-------------------!
	        	sp=0.
	        	su=g_const*alpha_l*(T_tent(i,j)-T_amb)*volume
	
	      		!----------------------------------!
				! Establish neighbour coefficients !
				!----------------------------------!
	        	aw(i,j)=dw+fw/2
	        	ae(i,j)=de-fe/2
	        	as(i,j)=ds+fs/2
	        	an(i,j)=dn-fn/2
	        
	
	      		!------------------------------!
				! Establish center coefficient !
				!------------------------------!
	        	ap_zero= rho*volume/dt
	        	ap_tilde=aw(i,j)+ae(i,j)+as(i,j)+an(i,j)+fe-fw+fn-fs-sp
	        	ap(i,j)=ap_zero+ap_tilde
	
	      		!-----------------------------------------------!
				! Calculates coefficient term for use in pc eq. !
				!-----------------------------------------------!
	        	d_v(i,j)=areas*relax(2)/ap(i,j)
	
	      		!-----------------------------!
				! Establish total source term !
				!-----------------------------!
	        	b(i,j)=(p(i,j-1)-p(i,j))*areas+ap_zero*v(i,j)+su
	
	      		!--------------------------------------------------------------!
				! Introducing relaxation and adjusting source term accordingly !
				!--------------------------------------------------------------!
	        	ap(i,j)=ap(i,j)/relax(1)
	        	b(i,j)=b(i,j)+(1-relax(1))*ap(i,j)*v_tent(i,j)
	      	end do
	   	end do
	
	end subroutine vcoeff
	
	
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!-------------------------------------------------------------------------!
	!-------------Calculates coefficients for the temp. equation--------------!
	!-------------------------------------------------------------------------!
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	subroutine tcoeff(ae,aw,an,as,ap,b)
		
		!----------------!
		! Initialization !
		!----------------!
	  	use declarations
	  	implicit none
	  	double precision, dimension(:,:), intent(out) :: ae,aw,an,as,ap,b
	  	integer :: i,j
	  	real :: fw,fe,fs,fn,dw,de,ds,dn,areaw,areae,areas,arean,sp,su
	  	real :: ap_tilde,ap_zero,h_tot
	    real :: r_w,r_e,dx_cell,dy_cell,volume
	    
		!-----------------------------!
		! Calculates convective terms !
		!-----------------------------!
	  	call conv()
		
		!-----------------------!
		! Iterate through array !
		!-----------------------!
	  	do i=2,npi-1
	    	do j=2,npj-1
				
				!-------------------------!
				! Find area of cell faces !
				!-------------------------!
!	      		areaw=y_v(j+1)-y_v(j)
!	      		areae=areaw
!	      		areas=x_u(i+1)-x_u(i)
!	      		arean=areas
	      		dx_cell=y_v(j+1)-y_v(j)
				dy_cell=x_u(i+1)-x_u(i)
				r_w=abs(xl/2-x_u(i))
				r_e=abs(xl/2-x_u(i+1))
				
				areaw=pi*r_w*dy_cell
	      		areae=pi*r_e*dy_cell
	      		areas=pi/2*abs(r_w**2-r_e**2)
	      		arean=areas
	      		volume=areas*dy
				
				!---------------------------!
				! Find convective mass flux !
				!---------------------------!
	      		fw=f_u(i,j)*cp*areaw
	      		fe=f_u(i+1,j)*cp*areae
	      		fs=f_v(i,j)*cp*areas
	      		fn=f_v(i,j+1)*cp*arean

				!---------------------------------------------!
				! Find diffusion conductance by harmonic mean !
				!---------------------------------------------!
!				dw=((k(i-1,j)*k(i,j))/(k(i-1,j)*    &
!	            	(x(i)-x_u(i))+k(i,j)*(x_u(i)-x(i-1))))*areaw
!	        	de=((k(i,j)*k(i+1,j))/(k(i,j)*       &
!	            	(x(i+1)-x_u(i+1))+k(i+1,j)*(x_u(i+1)-x(i))))*areae
!	        	ds=((k(i,j-1)*k(i,j))/(k(i,j-1)*      &
!	            	(y(j)-y_v(j))+k(i,j)*(y_v(j)-y(j-1))))*areas
!	        	dn=((k(i,j)*k(i,j+1))/(k(i,j)*         &
!	            	(y(j+1)-y_v(j+1))+k(i,j+1)*(y_v(j+1)-y(j))))*arean
				dw=k/(x(i)-x(i-1))*areaw		
				de=k/(x(i+1)-x(i))*areae
				ds=k/(y(j)-y(j-1))*areas
				dn=k/(y(j+1)-y(j))*arean
	
				!-------------------!
				! Find source terms !
				!-------------------!
				sp=0.
	      		su=0.
	      		! Check if next to wall. If true, add heat loss through wall
	        	if (i==2) then
	        		Ra=g_const/(nu_a*beta_a*T_amb)*(T(i,npj/2)-T_amb)*yl**3
	        		h_co=k_a/yl*( 0.68+0.67*Ra**0.25/	&
	        			( 1+(0.492/Pr_a)**(9/16) )**(4/9) )    		
	        		h_tot=1/(cup_width/k_cup + 1/h_co)
	        		
	          		su=su-h_tot*areaw*(T(i,j)-T_amb)
	        	elseif (i==npi-1) then
	        		Ra=g_const/(nu_a*beta_a*T_amb)*(T(i,npj/2)-T_amb)*yl**3
	        		h_co=k_a/yl*( 0.68+0.67*Ra**0.25/	&
	        			( 1+(0.492/Pr_a)**(9/16) )**(4/9) )
	        		h_tot=1/(cup_width/k_cup + 1/h_co)
	        		
	        		su=su-h_tot*areae*(T(i,j)-T_amb)
	        	end if
	        	if (j==2) then
	        		!Ra=g_const/(nu_a*beta_a*T(i,j))*(T(i,j)-T_amb)*areaw**3
	        		!h_co=k_a*0.54*Ra**0.25/arean
	        		h_co=10
	        		h_tot=1/(cup_width/k_cup + 1/h_co)
	        		
	          		su=su-h_tot*areas*(T(i,j)-T_amb)
	        	elseif (j==npj-1) then
	        		!Ra=g_const/(nu_a*beta_a*T(i,j))*(T(i,j)-T_amb)*areaw**3
	        		!h_co=k_a*0.27*Ra**0.25/areas
	        		h_co=45
	        		h_tot=1/(cup_width/k_cup + 1/h_co)
	        		
	        		su=su-h_tot*arean*(T(i,j)-T_amb)
	        	end if
	
	            
	     		!----------------------------------!
				! Establish neighbour coefficients !
				!----------------------------------!
	        	aw(i,j)=dw+fw/2
	        	ae(i,j)=de-fe/2
	        	as(i,j)=ds+fs/2
	        	an(i,j)=dn-fn/2            
	            
	            !------------------------------!
				! Establish center coefficient !
				!------------------------------!
	        	ap_zero= rho*cp*volume/dt
	        	ap_tilde=aw(i,j)+ae(i,j)+as(i,j)+an(i,j)+fe-fw+fn-fs-sp
	        	ap(i,j)=ap_zero+ap_tilde
	      	
				!-----------------------------!
				! Establish total source term !
				!-----------------------------!
	     	 	b(i,j)=ap_zero*T(i,j)+su
	    	end do
	  	end do   
	  	
	end subroutine tcoeff
	
	
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!-------------------------------------------------------------------------!
	!--------------------Calculates convective mass flux----------------------!
	!-------------------------------------------------------------------------!
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	subroutine conv()
		
		!----------------!
		! Initialization !
		!----------------!     
	  	use declarations
	  	implicit none
	  	integer :: i,j
		
		!---------------------------------------------!
		! Calculates conv. mass flux for whole domain !
		!---------------------------------------------!
	  	do i=2,npi
	    	do j=2,npj
	      		f_u(i,j)= (rho*(x(i)-x_u(i))+           &
	            	rho*(x_u(i)-x(i-1)) )*u_tent(i,j)/(x(i)-x(i-1))
	      		f_v(i,j)= (rho*(y(j)-y_v(j))+           &
	            	rho*(y_v(j)-y(j-1)))*v_tent(i,j)/(y(j)-y(j-1))
	    	end do
	 	end do   
	     
	end subroutine conv
	
	
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!-------------------------------------------------------------------------!
	!----------------Calculates coefficients for the pc equ.------------------!
	!-------------------------------------------------------------------------!
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	subroutine pccoeff(ae,aw,an,as,ap,b)
		
		!----------------!
		! Initialization !
		!----------------!
	  	use declarations
	 	implicit none
	  	double precision, dimension(:,:), intent(out) :: ae,aw,an,as,ap,b
	  	integer :: i,j
	  	real ::  areaw,areae,areas,arean,sp
	  	real :: r_w,r_e,dx_cell,dy_cell,volume
		
		!-----------------------------!
		! Calculates convective terms !
		!-----------------------------!
	  	call conv()
		
		!-----------------------!
		! Iterate through array !
		!-----------------------!
	  	do i=2,npi-1
	    	do j=2,npj-1
	    	
				!-------------------------!
				! Find area of cell faces !
				!-------------------------!
!	      		areaw=y_v(j+1)-y_v(j)
!	      		areae=areaw
!	      		areas=x_u(i+1)-x_u(i)
!	      		arean=areas
	      		dx_cell=y_v(j+1)-y_v(j)
				dy_cell=x_u(i+1)-x_u(i)
				r_w=abs(xl/2-x_u(i))
				r_e=abs(xl/2-x_u(i+1))
				
				areaw=pi*r_w*dy_cell
	      		areae=pi*r_e*dy_cell
	      		areas=pi/2*abs(r_w**2-r_e**2)
	      		arean=areas
	      		volume=areas*dy
				
				!-------------------!
				! Find source terms !
				!-------------------!
	      		b(i,j)=f_u(i,j)*areaw-f_u(i+1,j)*areae+   &
	         		f_v(i,j)*areas-f_v(i,j+1)*arean
	      		sp=0.   
				
				!----------------------------------!
				! Establish neighbour coefficients !
				!----------------------------------!              
	      		ae(i,j)=rho*d_u(i+1,j)*areae
	      		aw(i,j)=rho*d_u(i,j)*areaw
	      		an(i,j)=rho*d_v(i,j+1)*arean
	      		as(i,j)=rho*d_v(i,j)*areas
				
				!------------------------------!
				! Establish center coefficient !
				!------------------------------!
	      		ap(i,j)=ae(i,j)+aw(i,j)+an(i,j)+as(i,j)-sp
	      		
	      		!----------------------------------------------------------!
				! Set old pc-values to zero to avoid influence on new ones !
				!----------------------------------------------------------!
	      		pc(i,j)=0.	
	    	end do
	  	end do
	end subroutine pccoeff
	
	
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!-------------------------------------------------------------------------!
	!----------------Calculates the 2-norm error of an array------------------!
	!-------------------------------------------------------------------------!
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	subroutine TwoNormError(array,i)
	 	
	 	!----------------!
		! Initialization !
		!----------------!
		use declarations
	  	implicit none
	  	double precision, dimension(:,:), intent(in) :: array	
	  	integer :: i
	  
	  	!-----------------!
		! Calculate error !
		!-----------------!
	  	er(i) = ( dx*dy*sum( (array)**2) )**0.5
	  
	end subroutine TwoNormError

end module sub


!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
!----------------------------------------------------------------------------!
! SOLVES TRANSIENT CONVECTION-DIFFUSION PROBLEMS USING THE SIMPLE ALGORITHM  !
!----------------------------------------------------------------------------!
! Calculates velocities, pressure and temperature development in time on a   !
! 2D cartesian grid.                                                         !
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
 program SIMPLE
		
	!----------------!
	! Initialization !
	!----------------!
	use declarations
	use sub
	implicit none
	double precision,allocatable,dimension(:,:) :: ae,aw,an,as,ap,b
	integer :: steps,it,jt,n
	integer :: txtclock, coeffclock, solveclock
  	real    :: txttime, coefftime, solvetime
  	
  	!-------------------------------------------!
	! Initialize all parameters and create grid !
	!-------------------------------------------!
 	call init()
 	call grid()
 	
 	!----------------------------------------------------------!
	! Specify BCs that remain untouched through the iterations !
	!----------------------------------------------------------!    
  	call fixedbound()
 
	!---------------------------!
	! Allocate remaining arrays !
	!---------------------------!
	allocate(ae(npi,npj),aw(npi,npj),an(npi,npj),as(npi,npj))
  	allocate(ap(npi,npj),b(npi,npj))
  	allocate(f_u(npi,npj),f_v(npi,npj))

	!----------------------------------!
	! Set node for convergence history !
	!----------------------------------!
  	it=npi/2
  	jt=npj/2
  	write(*,*) 'Node for convergence history:',it,jt

	!--------------------------------------------------!
	! Set tentative values equal to initialized values !
	!--------------------------------------------------!
	u_tent=u
 	v_tent=v
 	T_tent=T
 	
 	!---------------!
	! March in time !
	!---------------!
	do n = 1,nsteps
		call tick(txtclock)	 
		!----------------------------------------------------------------!
		! Iterate 'last' nr. of steps or until convergence criterion met !
		!----------------------------------------------------------------!		    
	    do iter=1,last
	      
	      	!--------------------------------------!
			! Solve u-equation for tentative value !
			!--------------------------------------!
			call tick(coeffclock)
	      	call ucoeff(ae,aw,an,as,ap,b)
	      	coefftime=tock(coeffclock)
	      	!print *, 'coeff:', coefftime
	      	call tick(solveclock)
	    	do steps=1,1
				call solve(u_tent,b,ae,aw,an,as,ap,2,npi,1,npj)
			end do
			solvetime=tock(solveclock)
			!print *, 'solve:', solvetime
			
			!--------------------------------------!
			! Solve v-equation for tentative value !
			!--------------------------------------!
			call vcoeff(ae,aw,an,as,ap,b)
			do steps=1,1
				call solve(v_tent,b,ae,aw,an,as,ap,1,npi,2,npj)
			end do
			
			!-------------------!
			! Solve pc-equation !
			!-------------------!
			call pccoeff(ae,aw,an,as,ap,b)
			do steps=1,15
				call solve(pc,b,ae,aw,an,as,ap,1,npi,1,npj)
			end do
			
			!---------------------------------!
			! Correct velocities and pressure !
			!---------------------------------!
			call velcorr()	
			
			!--------------------------------------------------------------!
			! Calculate 2-norm error of conserved mass, i.e. b from pc-eq. !
			!--------------------------------------------------------------!	
			call TwoNormError(b,iter)
			
			!---------------------!
			! Solve temp-equation !
			!---------------------!
			call tcoeff(ae,aw,an,as,ap,b)  
			do steps=1,1
				call solve(T_tent,b,ae,aw,an,as,ap,1,npi,1,npj)
			end do
			
			!-------------------------!
			! Specify boundary values !
			!-------------------------!
			call bound()			
			
!			write(*,*) '(er(iter))/(er(1): ',(er(iter))/er(1)
			!-----------------------!
			! Check for convergence !
			!-----------------------!
			if((er(iter))/(er(1))<=0.0001)then
			
				!-------------------------------------!
				! Convergence met -> break iterations !
				!-------------------------------------!
				write(*,*) 'needed iterations:',iter
				last=iter
				exit
				
							
			else if(iter==last) then
				
				!------------------------------------!
				! Convergence not met -> print error !
				!------------------------------------!
				write(*,*) 'No convergence, er: ',er(iter)/er(1)
				
			end if	
		end do
		
		!---------------!
		! Update values !
		!---------------!
		u=u_tent
		v=v_tent
		T=T_tent

		!---------------------------------------!
		! Print results for every 10th timestep !
		!---------------------------------------!
   		if (mod(n,1) == 0) then
      		write (*,'(f7.2,5g15.5)')  n*dt, &
             	u(it,jt),v(it,jt),p(it,jt),pc(it,jt),T(it,jt)
            
        	call print(real(dt*n)) 
        	
    	end if
    	txttime=tock(txtclock)
    	print *, 'System time for timestep = ', txttime  
    		    	      
	end do
     
end program SIMPLE

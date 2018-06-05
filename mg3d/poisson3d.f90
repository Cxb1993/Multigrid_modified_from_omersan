!   Modified from  https://github.com/omersan/5.06.Multigrid2D
program poisson3d
    implicit none
    integer::i,j,k
    integer::nx,ny,nz
    integer::isolver
    real(kind=8) ::x0,xL
    real(kind=8) ::y0,yL
    real(kind=8) ::z0,zL
    real(kind=8) ::dx,dy,dz
    real(kind=8) ::tol,rms
    real(kind=8),dimension(:,:,:),allocatable ::u,f,ue,e
    real(kind=8),dimension(:),allocatable ::x,y,z



    !Domain
    x0 =-1.0d0 !left
    xL = 1.0d0 !right

    y0 =-1.0d0 !bottom
    yL = 1.0d0 !up

    z0 =-1.0d0 !bottom
    zL = 1.0d0 !up

    !number of points
    nx = 64 !number of grid points in x (i.e., should be power of 2)
    ny = nx  !number of grid points in y
    nz = nx 

    !grid spacing (spatial)
    dx = (xL-x0)/dfloat(nx)
    dy = (yL-y0)/dfloat(ny)
    dz = (zL-z0)/dfloat(nz)

    !spatial coordinates 
    allocate(x(0:nx))
    do i=0,nx
        x(i) = x0 + dfloat(i)*dx
    end do

    allocate(y(0:ny))
    do j=0,ny
        y(j) = y0 + dfloat(j)*dy
    end do


    allocate(z(0:nz))
    do k=0,nz
        z(k) = z0 + dfloat(k)*dz
    end do


    !Tolerance
    tol= 1.0d-7


    allocate( u(0:nx,0:ny,0:nz))
    allocate( f(0:nx,0:ny,0:nz))
    allocate( e(0:nx,0:ny,0:nz))
    allocate(ue(0:nx,0:ny,0:nz))

!---------------------------------------------!
!   Exact solution (test case from Moin's textbook):
!---------------------------------------------!
    do k=0,nz
    do j=0,ny
    do i=0,nx
         f(i,j,k) = 2.0d0 * ( (x(i)*x(i)-1.0d0)*(y(j)*y(j)-1.0d0)	&
			+   (y(j)*y(j)-1.0d0)*(z(k)*z(k)-1.0d0)		&
			+   (z(k)*z(k)-1.0d0)*(x(i)*x(i)-1.0d0)  )


        ue(i,j,k)= (x(i)*x(i)-1.0d0)*(y(j)*y(j)-1.0d0)*(z(k)*z(k)-1.0d0)
    end do
    end do
    end do

    !Initial guess Numerical solution:
    do i=0,nx
    do j=0,ny
    do k=0,nz
        u(i,j,k)=0.0d0
    end do
    end do
    end do

    !Boundary conditions has to satisfy exact solution

        u(0 ,:,:)  = ue(0 ,:,:)		
        u(nx,:,:)  = ue(nx,:,:)	

        u(:,0 ,:)  = ue(:,0 ,:)	
        u(:,ny,:)  = ue(:,ny,:)					  					  	

        u(:,:,0 )  = ue(:,:,0 )		
        u(:,:,nz)  = ue(:,:,nz)						  	



    open(19,file='output.txt')

!----------------------!
!   Solver:
!----------------------!
    isolver = 0
    if (isolver.eq.0) then !Gauss-Seidel scheme
	 call gauss_seidel(nx,ny,nz,dx,dy,dz,f,u,tol)
    else
	    call multigrid(nx,ny,nz,dx,dy,dz,f,u,tol)
    end if


!----------------------!
!   Error analysis:
!----------------------!
    do i=0,nx
    do j=0,ny
    do k=0,nz
        e(i,j,k) = dabs(u(i,j,k)-ue(i,j,k))
    end do 
    end do
    end do

    !L-2 Norm:
    call l2norm(nx,ny,nz,e,rms)

    write(*,*)"L2-norm =",rms
    write(19,*)"L2-norm =",rms

    !maximum norm
    write(*,*)"Max-norm =",maxval(e)
    write(19,*)"Max-norm =",maxval(e)
    close(19)

    !Plot field
    open(10,file='field.plt')
    write(10,*) 'variables ="x","y","z","f","u","ue"'
    write(10,*)'zone f=point i=',nx+1,',j=',ny+1,',k=',nz+1
    do k=0,nz
    do j=0,ny
    do i=0,nx
        write(10,*) x(i),y(j),z(k),f(i,j,k),u(i,j,k),ue(i,j,k)
    end do
    end do
    end do
    close(10)


    end program poisson3d    

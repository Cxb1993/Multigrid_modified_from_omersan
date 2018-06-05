!     Modified from  https://github.com/omersan/5.06.Multigrid2D
!     V-cycle Multigrid method for solving 2D Poisson equation
!     d2u/dx2 + d2u/dy2 = f(x,y)
!     Drichlet b.c.

program poisson2d
    implicit none
    integer::i,j,nx,ny,isolver
    real(kind=8),dimension(:,:),allocatable ::u,f,ue,e
    real(kind=8),dimension(:),allocatable ::x,y
    real(kind=8) ::dx,dy,tol,rms,x0,xL,y0,yL

    !Domain
    x0 =-1.0d0 !left
    xL = 1.0d0 !right

    y0 =-1.0d0 !bottom
    yL = 1.0d0 !up

    !number of points
    nx = 128 !number of grid points in x (i.e., should be power of 2)
    ny = nx  !number of grid points in y

    !grid spacing (spatial)
    dx = (xL-x0)/dfloat(nx)
    dy = (yL-y0)/dfloat(ny)

    !spatial coordinates 
    allocate(x(0:nx))
    do i=0,nx
        x(i) = x0 + dfloat(i)*dx
    end do

    allocate(y(0:ny))
    do j=0,ny
        y(j) = y0 + dfloat(j)*dy
    end do





    !Tolerance
    tol= 1.0d-7


    allocate(u(0:nx,0:ny))
    allocate(f(0:nx,0:ny))
    allocate(e(0:nx,0:ny))
    allocate(ue(0:nx,0:ny))

!---------------------------------------------!
!   Exact solution (test case from Moin's textbook):
!---------------------------------------------!
    do j=0,ny
    do i=0,nx
        f(i,j) =-2.0d0*(2.0d0-x(i)*x(i)-y(j)*y(j))
        ue(i,j)= (x(i)*x(i)-1.0d0)*(y(j)*y(j)-1.0d0)
    end do
    end do


    !Initial guess Numerical solution:
    do i=0,nx
    do j=0,ny
        u(i,j)=0.0d0
    end do
    end do

    !Boundary conditions has to satisfy exact solution
        u(:,0)  = ue(:,0)	
        u(:,ny) = ue(:,ny)					  					  	


        u(0,:)  = ue(0,:)		
        u(nx,:) = ue(nx,:)						  	


    open(19,file='output.txt')

!----------------------!
!   Solver:
!----------------------!
    isolver = 1
    if (isolver.eq.0) then !Gauss-Seidel scheme
	    call gauss_seidel(nx,ny,dx,dy,f,u,tol)
    else
	    call multigrid(nx,ny,dx,dy,f,u,tol)
    end if


!----------------------!
!   Error analysis:
!----------------------!
    do i=0,nx
    do j=0,ny
        e(i,j) = dabs(u(i,j)-ue(i,j))
    end do 
    end do


    !L-2 Norm:
    call l2norm(nx,ny,e,rms)

    write(*,*)"L2-norm =",rms
    write(19,*)"L2-norm =",rms

    !maximum norm
    write(*,*)"Max-norm =",maxval(e)
    write(19,*)"Max-norm =",maxval(e)
    close(19)

    !Plot field
    open(10,file='field.plt')
    write(10,*) 'variables ="x","y","f","u","ue"'
    write(10,*)'zone f=point i=',nx+1,',j=',ny+1
    do j=0,ny
    do i=0,nx
        write(10,*) x(i),y(j),f(i,j),u(i,j),ue(i,j)
    end do
    end do
    close(10)


    end program poisson2d    

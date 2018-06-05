!---------------------------------------------------------------------------!
!Residual formula for Poisson equation
!Works for Drichlet boundary conditions (Boundary points never updated)
!---------------------------------------------------------------------------!
subroutine residual(nx,ny,dx,dy,f,u,r)

    implicit none

    integer::nx,ny
    real(kind=8) ::dx,dy
    real(kind=8), dimension(0:nx,0:ny)::u,f,r
    integer::i,j

    do i=1,nx-1
    do j=1,ny-1

	r(i,j) = f(i,j) - (u(i+1,j) - 2.0d0*u(i,j) + u(i-1,j))/(dx*dx) &
		        - (u(i,j+1) - 2.0d0*u(i,j) + u(i,j-1))/(dy*dy) 

    end do
    end do

    !Boundary conditions for residuals
    do i=0,nx
	r(i,0)  = 0.0d0	
	r(i,ny) = 0.0d0					  					  	
    end do

    do j=0,ny
	r(0,j)  = 0.0d0		
	r(nx,j) = 0.0d0						  	
    end do

    return

end subroutine residual

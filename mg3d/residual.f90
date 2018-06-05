!---------------------------------------------------------------------------!
!Residual formula for Poisson equation
!Works for Drichlet boundary conditions (Boundary points never updated)
!---------------------------------------------------------------------------!
subroutine residual(nx,ny,nz,dx,dy,dz,f,u,r)

    implicit none

    integer	::nx,ny,nz
    real(kind=8)::dx,dy,dz
    real(kind=8), dimension(0:nx,0:ny,0:nz)::u,f,r
    integer::i,j,k

    do i=1,nx-1
    do j=1,ny-1
    do k=1,nz-1

	r(i,j,k) = f(i,j,k) - (u(i+1,j,k) - 2.0d0*u(i,j,k) + u(i-1,j,k))/(dx*dx)	&
		            - (u(i,j+1,k) - 2.0d0*u(i,j,k) + u(i,j-1,k))/(dy*dy)	&
		            - (u(i,j,k+1) - 2.0d0*u(i,j,k) + u(i,j,k-1))/(dy*dy)  

    end do
    end do
    end do

    !Boundary conditions for residuals

	r(0 ,:,:) = 0.0d0		
	r(nx,:,:) = 0.0d0						  	

	r(:,0 ,:) = 0.0d0	
	r(:,ny,:) = 0.0d0					  					  	

	r(:,:,0 ) = 0.0d0		
	r(:,:,nz) = 0.0d0						  	

    return

end subroutine residual

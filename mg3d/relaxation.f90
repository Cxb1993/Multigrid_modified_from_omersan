subroutine relaxation_sor(nx,ny,nz,dx,dy,dz,rhs,x_soln)

    implicit none

    integer					::	nx,ny,nz
    real(kind=8) 				::	omega=1.1
    real(kind=8) 				::	dx,dy,dz
    real(kind=8), dimension(0:nx,0:ny,0:nz)	::	x_soln,rhs
    real(kind=8) 				::	ae,aw,an,as,at,ab,ap
    integer					::	i,j,k

!keep it A-matrix format
 
     ae 	= 1./(dx*dx)
     aw 	= 1./(dx*dx)
     an 	= 1./(dy*dy) 
     as 	= 1./(dy*dy)
     at 	= 1./(dz*dz)
     ab 	= 1./(dz*dz)
     ap 	= -1.*(ae + aw + an + as + at + ab)

 
    do i=1,nx-1
    do j=1,ny-1
    do k=1,nz-1
!------------------SOR-type------------------------------------------------
  
	x_soln(i,j,k) = (omega*rhs(i,j,k)						&
			  - omega*( ae*x_soln(i+1,j  ,k  ) + aw*x_soln(i-1,j  ,k  ) 	&
				  + an*x_soln(i  ,j+1,k  ) + as*x_soln(i  ,j-1,k  ) 	&
				  + at*x_soln(i  ,j  ,k+1) + ab*x_soln(i  ,j  ,k-1) ) 	&
			  - (omega-1)*ap*x_soln(i,j,k) )/ap

!-------------------------------------------------------------------------
    end do
    end do
    end do

    return

end subroutine relaxation_sor

subroutine relaxation_gs(nx,ny,nz,dx,dy,dz,rhs,x_soln)

    implicit none

    integer					::	nx,ny,nz
    real(kind=8) 				::	dx,dy,dz
    real(kind=8), dimension(0:nx,0:ny,0:nz)	::	x_soln,rhs
    real(kind=8) 				::	ae,aw,an,as,at,ab,ap
    integer					::	i,j,k


!keep it A-matrix format

     ae 	= 1./(dx*dx)
     aw 	= 1./(dx*dx)
     an 	= 1./(dy*dy) 
     as 	= 1./(dy*dy)
     at 	= 1./(dz*dz)
     ab 	= 1./(dz*dz)
     ap 	= -1.*(ae + aw + an + as + at + ab)

  
    do i=1,nx-1
    do j=1,ny-1
    do k=1,nz-1
!------------------GS-type------------------------------------------------
        x_soln(i,j,k) =   (rhs(i,j,k) 							&
				-(  ae*x_soln(i+1,j  ,k  ) + aw*x_soln(i-1,j  ,k  ) 	&
				  + an*x_soln(i  ,j+1,k  ) + as*x_soln(i  ,j-1,k  ) 	&
				  + at*x_soln(i  ,j  ,k+1) + ab*x_soln(i  ,j  ,k-1) ) 	&
		  	)/ap
!-------------------------------------------------------------------------
    end do
    end do
    end do

    return

end subroutine relaxation_gs


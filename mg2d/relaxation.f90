subroutine relaxation_sor(nx,ny,dx,dy,rhs,x_soln)

    implicit none

    integer					::	nx,ny
    real(kind=8) 				::	omega=1.1
    real(kind=8) 				::	dx,dy
    real(kind=8), dimension(0:nx,0:ny)		::	x_soln,rhs
    real(kind=8) 				::	a_e,a_w,a_n,a_s,a_p
    integer					::	i,j

!keep it A-matrix format
 
     a_e 	= 1./(dx*dx)
     a_w 	= 1./(dx*dx)
     a_n 	= 1./(dy*dy) 
     a_s 	= 1./(dy*dy)
     a_p 	= -1.*(a_e + a_w + a_n + a_s)

  
    do i=1,nx-1
    do j=1,ny-1
!------------------SOR-type------------------------------------------------
  
	x_soln(i,j) = (omega*rhs(i,j)			&
			  - omega*( a_e*x_soln(i+1,j) 	&
				  + a_n*x_soln(i,j+1) 	&
				  + a_w*x_soln(i-1,j) 	&
				  + a_s*x_soln(i,j-1) ) &
			  - (omega-1)*a_p*x_soln(i,j)	&
			)/a_p

!-------------------------------------------------------------------------
    end do
    end do

    return

end subroutine relaxation_sor

subroutine relaxation_gs(nx,ny,dx,dy,rhs,x_soln)

    implicit none

    integer					::	nx,ny
    real(kind=8) 				::	dx,dy
    real(kind=8), dimension(0:nx,0:ny)		::	x_soln,rhs
    real(kind=8) 				::	a_e,a_w,a_n,a_s,a_p
    integer					::	i,j


!keep it A-matrix format

     a_e 	= 1./(dx*dx)
     a_w 	= 1./(dx*dx)
     a_n 	= 1./(dy*dy) 
     a_s 	= 1./(dy*dy)
     a_p 	= -1.*(a_e + a_w + a_n + a_s)


  
    do i=1,nx-1
    do j=1,ny-1
!------------------GS-type------------------------------------------------
        x_soln(i,j) =   (rhs(i,j) 			&
				- a_e*x_soln(i+1,j)	&
				- a_w*x_soln(i-1,j)	&
				- a_n*x_soln(i,j+1)	&
				- a_s*x_soln(i,j-1)	&
		  	)/a_p
!-------------------------------------------------------------------------
    end do
    end do

    return

end subroutine relaxation_gs


!---------------------------------------------------------------------------!
!Restriction operators
!---------------------------------------------------------------------------!


subroutine restriction(nxc,nyc,nxf,nyf,val_c,val_f)

    implicit none

    integer					:: nxc,nyc
    integer					:: nxf,nyf
    real(kind=8), dimension(0:nxc,0:nyc)	:: val_c	!on coarse grid
    real(kind=8), dimension(0:nxf,0:nyf)	:: val_f	!on fine grid
    integer					:: i,j
    integer					:: ireo

    ireo = 3

    if (ireo.eq.1) then !simply injection

	do i=1,nxf-1
	do j=1,nyf-1
	    val_f(i,j) = val_c(2*i,2*j) 							  	
	end do
	end do

    else if (ireo.eq.2) then !half-weight

	do i=1,nxf-1
	do j=1,nyf-1
		val_f(i,j) = 1.0d0/8.0d0*( 4.0d0*val_c(2*i,2*j) &
				+ 1.0d0*(val_c(2*i+1,2*j)	&
					+val_c(2*i-1,2*j)	&
					+val_c(2*i,2*j+1)	&
					+val_c(2*i,2*j-1)) )							  	
	end do
	end do


    else !full-weight (trapezoidal)

	do i=1,nxf-1
	do j=1,nyf-1
	    val_f(i,j) = (1.0d0/16.0d0)*( 4.0d0*val_c(2*i,2*j) &
	     		+ 2.0d0*(val_c(2*i+1,2*j)		&
				+val_c(2*i-1,2*j)		&
				+val_c(2*i,2*j+1)		&
				+val_c(2*i,2*j-1)) 		&
	     		+ 1.0d0*(val_c(2*i+1,2*j+1)		&
				+val_c(2*i-1,2*j-1)		&
				+val_c(2*i-1,2*j+1)		&
				+val_c(2*i+1,2*j-1)))							  	
	end do
	end do

    end if


    !update boundaries
    do i=0,nxf
	val_f(i,0)   = val_c(2*i,0) 	
	val_f(i,nyf) = val_c(2*i,nyc) 					  					  	
    end do

    do j=0,nyf
	val_f(0,j)   = val_c(0,2*j)		
	val_f(nxf,j) = val_c(nxc,2*j)						  	
    end do

  
    return

end subroutine restriction

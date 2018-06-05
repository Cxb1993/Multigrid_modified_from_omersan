!---------------------------------------------------------------------------!
!Restriction operators
!---------------------------------------------------------------------------!


subroutine restriction(nxc,nyc,nzc,nxf,nyf,nzf,val_c,val_f)

    implicit none

    integer					:: nxc,nyc,nzc
    integer					:: nxf,nyf,nzf
    real(kind=8), dimension(0:nxc,0:nyc,0:nzc)	:: val_c	!on coarse grid
    real(kind=8), dimension(0:nxf,0:nyf,0:nzf)	:: val_f	!on fine grid
    integer					:: i,j,k, ic,jc,kc
    integer					:: ireo

    ireo = 3

    if (ireo.eq.1) then !simply injection

	do i=1,nxf-1
	do j=1,nyf-1
	do k=1,nzf-1
            
		ic=2*i
		jc=2*j
		kc=2*k
		val_f(i,j,k) = val_c(ic,jc,kc) 							  	
	end do
	end do
	end do

    else !full-weight (trapezoidal)

	do i=1,nxf-1
	do j=1,nyf-1
	do k=1,nzf-1
            
        ic=2*i
        jc=2*j
        kc=2*k
	val_f(i,j,k) = (1.0d0/84.0d0)*( 	 12.0d0* 	 val_c(ic,jc,kc) 		&
					+ 4.0d0*(  	 val_c(ic+1,jc,kc)		&
							+val_c(ic-1,jc,kc)		&
		  					+val_c(ic,jc+1,kc)		&
							+val_c(ic,jc-1,kc)		&
		  					+val_c(ic,jc,kc+1)		&
							+val_c(ic,jc,kc-1) ) 		&
					+ 3.0d0*(  	 val_c(ic+1,jc+1,kc)		&
							+val_c(ic+1,jc-1,kc)		&
							+val_c(ic-1,jc+1,kc)		&
							+val_c(ic-1,jc-1,kc)		&
		  					+val_c(ic,jc+1,kc+1)		&
							+val_c(ic,jc+1,kc-1)		&
							+val_c(ic,jc-1,kc+1)		&
							+val_c(ic,jc-1,kc-1)		&
   		  					+val_c(ic+1,jc,kc+1)		&
							+val_c(ic+1,jc,kc-1)		&
							+val_c(ic-1,jc,kc+1)		&
							+val_c(ic-1,jc,kc-1) 	)	&
					+ 2.0d0*(	 val_c(ic+1,jc+1,kc+1)		&
							+val_c(ic+1,jc-1,kc-1)		&
							+val_c(ic+1,jc-1,kc+1)		&
							+val_c(ic+1,jc+1,kc-1)		&
							+val_c(ic-1,jc+1,kc+1)		&
							+val_c(ic-1,jc+1,kc-1)		&
							+val_c(ic-1,jc-1,kc+1)		&
							+val_c(ic-1,jc-1,kc-1)	))							  	
	end do
	end do
	end do

    end if


    !update boundaries
    do j=0,nyf
    do k=0,nzf
	jc=2*j
	kc=2*k
	val_f(0  ,j,k) = val_c(0  ,jc,kc)		
	val_f(nxf,j,k) = val_c(nxc,jc,kc)						  	
    end do
    end do

    do k=0,nzf
    do i=0,nxf
	kc=2*k
	ic=2*i
	val_f(i,0  ,k) = val_c(ic,0  ,kc) 	
	val_f(i,nyf,k) = val_c(ic,nyc,kc) 					  					  	
    end do
    end do


    do i=0,nxf
    do j=0,nyf
	ic=2*i
	jc=2*j
	val_f(i,j,0  ) = val_c(ic,jc,0  )		
	val_f(i,j,nzf) = val_c(ic,jc,nzc)						  	
    end do
    end do
  
    return

end subroutine restriction

!---------------------------------------------------------------------------!
!Prolongation operator
!bilinear interpolation
!---------------------------------------------------------------------------!

subroutine prolongation(nxf,nyf,nzf,nxc,nyc,nzc,val_f,val_c)

    implicit none

    integer					:: nxc,nyc,nzc
    integer					:: nxf,nyf,nzf
    real(kind=8), dimension(0:nxc,0:nyc,0:nzc)	:: val_c	!on coarse grid
    real(kind=8), dimension(0:nxf,0:nyf,0:nzf)	:: val_f	!on fine grid
    integer					:: i,j,k
    integer					:: ic,jc,kc


    do i=0,nxf-1
    do j=0,nyf-1
    do k=0,nzf-1

        ic=2*i
        jc=2*j
        kc=2*k

	val_c(ic  ,jc  ,kc  )	= val_f(i,j,k)

	val_c(ic+1,jc  ,kc  )	= (1.0d0/2.0d0)*(val_f(i,j,k)+val_f(i+1,j,k))
        val_c(ic  ,jc+1,kc  )	= (1.0d0/2.0d0)*(val_f(i,j,k)+val_f(i,j+1,k))
        val_c(ic  ,jc  ,kc+1)	= (1.0d0/2.0d0)*(val_f(i,j,k)+val_f(i,j,k+1))

	val_c(ic+1,jc+1,kc  )	= (1.0d0/4.0d0)*(val_f(i,j,k)+val_f(i+1,j,k)+val_f(i,j+1,k)+val_f(i+1,j+1,k))
	val_c(ic  ,jc+1,kc+1)	= (1.0d0/4.0d0)*(val_f(i,j,k)+val_f(i,j+1,k)+val_f(i,j,k+1)+val_f(i,j+1,k+1))
	val_c(ic+1,jc  ,kc+1)	= (1.0d0/4.0d0)*(val_f(i,j,k)+val_f(i,j,k+1)+val_f(i+1,j,k)+val_f(i+1,j,k+1))

	val_c(ic+1,jc  ,kc+1)	= (1.0d0/8.0d0)*(val_f(i,j,k)							&
						+ val_f(i+1,j,k)+val_f(i,j+1,k)+val_f(i,j,k+1)			&
						+ val_f(i+1,j+1,k)+val_f(i,j+1,k+1)+val_f(i+1,j,k+1)		&
						+ val_f(i+1,j+1,k+1) )

    end do
    end do
    end do

    do j=0,nyf
    do k=0,nzf
	val_c(nxc,2*j,2*k)    = val_f(nxf,j,k)
    end do
    end do

    do k=0,nzf
    do i=0,nxf
	val_c(2*i,nyc,2*k)    = val_f(i,nyf,k)
    end do
    end do

    do i=0,nxf
    do j=0,nyf
	val_c(2*i,2*j,nzc)    = val_f(i,j,nzf)
    end do
    end do


    return

end subroutine prolongation


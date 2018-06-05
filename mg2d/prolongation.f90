!---------------------------------------------------------------------------!
!Prolongation operator
!bilinear interpolation
!---------------------------------------------------------------------------!

subroutine prolongation(nxf,nyf,nxc,nyc,val_f,val_c)

    implicit none

    integer					:: nxc,nyc
    integer					:: nxf,nyf
    real(kind=8), dimension(0:nxc,0:nyc)	:: val_c	!on coarse grid
    real(kind=8), dimension(0:nxf,0:nyf)	:: val_f	!on fine grid
    integer					:: i,j


    do i=0,nxf-1
    do j=0,nyf-1
	val_c(2*i,2*j)    = val_f(i,j)
	val_c(2*i+1,2*j)  = 1.0d0/2.0d0*(val_f(i,j)+val_f(i+1,j))
        val_c(2*i,2*j+1)  = 1.0d0/2.0d0*(val_f(i,j)+val_f(i,j+1))
	val_c(2*i+1,2*j+1)= 1.0d0/4.0d0*(val_f(i,j)+val_f(i,j+1)+val_f(i+1,j)+val_f(i+1,j+1))
    end do
    end do

    do j=0,nyf
	val_c(nxc,2*j)    = val_f(nxf,j)
    end do

    do i=0,nxf
	val_c(2*i,nyc)    = val_f(i,nyf)
    end do

    return

end subroutine prolongation


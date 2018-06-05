!---------------------------------------------------------------------------!
!Compute L2-norm for an array
!---------------------------------------------------------------------------!
subroutine l2norm(nx,ny,nz,r,rms)

    implicit none

    integer	::nx,ny,nz
    real(kind=8), dimension(0:nx,0:ny,0:nz)::r
    integer	 ::i,j,k
    real(kind=8) ::rms,tmp

    tmp=0.0d0

    do i=1,nx-1
    do j=1,ny-1
    do k=1,nz-1
        tmp = tmp + r(i,j,k)*r(i,j,k)
    end do 
    end do
    end do

    rms= dsqrt(tmp/dfloat((nx-1)*(ny-1)*(nz-1)))

    return

end subroutine l2norm


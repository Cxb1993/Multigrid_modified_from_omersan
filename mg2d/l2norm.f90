!---------------------------------------------------------------------------!
!Compute L2-norm for an array
!---------------------------------------------------------------------------!
subroutine l2norm(nx,ny,r,rms)

    implicit none

    integer::Nx,Ny
    real(kind=8), dimension(0:nx,0:ny)::r
    integer::i,j
    real(kind=8) ::rms

    rms=0.0d0

    do i=1,nx-1
    do j=1,ny-1
        rms = rms + r(i,j)*r(i,j)
    end do 
    end do

    rms= dsqrt(rms/dfloat((nx-1)*(ny-1)))

    return

end subroutine l2norm


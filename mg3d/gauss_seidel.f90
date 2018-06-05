!---------------------------------------------------------------------------!
!Gauss Seidel scheme (1 level)
!---------------------------------------------------------------------------!
subroutine gauss_seidel(nx,ny,nz,dx,dy,dz,f,u,tol)

    implicit none

    integer		:: nx,ny,nz
    integer		:: k,ke,wl,n_iter
    real(kind=8) 	:: dx,dy,dz
    real(kind=8) 	:: tol
    real(kind=8) 	:: rms0,rms
    real(kind=8),dimension(0:nx,0:ny,0:nz)	:: u,f
    real(kind=8),dimension(:,:,:),allocatable	:: r



    n_iter = 100000  !maximum number of iteration

    allocate(r(0:nx,0:ny,0:nz))

    !Compute initial resitual and l2 norm:
    call residual(nx,ny,nz,dx,dy,dz,f,u,r)
    call l2norm(nx,ny,nz,r,rms0)

    open(66,file='residual.plt')
    write(66,*) 'variables ="k","rms","rms/rms0"'

    do k=1,n_iter

	    call relaxation_gs(nx,ny,nz,dx,dy,dz,f,u)

	    call residual(nx,ny,nz,dx,dy,dz,f,u,r)
    
	    ! Check for convergence on smallest grid	
	    call l2norm(nx,ny,nz,r,rms)

	    if (rms/rms0.le.tol) goto 10

            ! Write residual history
	    write(66,*) k,rms,rms/rms0

   	    write(*,*) k,rms,rms/rms0
     
    end do

10  continue
    
    close(66)

    deallocate(r)

    ke=k

    !work load (total number of operations)
    wl = ke*(nx*ny)

    write(19,*)"outer number of iteration = ",ke
    write(19,*)"normalized workload       = ",dfloat(wl)/dfloat(nx*ny)
    write(*,*)"outer number of iteration = ",ke
    write(*,*)"normalized workload       = ",dfloat(wl)/dfloat(nx*ny)

    return
 
end subroutine gauss_seidel  


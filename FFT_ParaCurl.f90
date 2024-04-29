! Program to find the Curl of a 3 dimensional scalar function by parallelly computing the FFT

module nodeinfo
    implicit none
    integer :: rank
    integer :: nprocs
end module nodeinfo

program fft_para_curl

    use nodeinfo
    use mpi

    implicit none

    include '/usr/local/include/fftw3.f'

    integer :: ierr
    integer, dimension(MPI_STATUS_SIZE) :: status1

    integer, parameter :: n = 250
    integer, parameter :: nx = n, ny = n, nz = n
    complex, parameter :: Z = (0,1)

    double complex, allocatable :: fhat1x(:), fhat1y(:), fhat1z(:), dfhat1x(:), dfhat1y(:), dfhat1z(:)

    real(8), allocatable :: actual_curl(:,:,:,:), curl(:,:,:,:), func(:,:,:,:)
    real(8), allocatable :: curltempx(:,:,:,:), curltempy(:,:,:,:), curltempz(:,:,:,:), curltemp(:,:,:,:)
    real(8), allocatable :: kappax(:), kappay(:), kappaz(:)
    real(8), allocatable :: der1x(:), der1y(:), der1z(:)
    real(8), allocatable :: derx(:,:,:,:), dery(:,:,:,:), derz(:,:,:,:), der(:,:,:,:)

    integer(8) :: forward1x, forward1y, forward1z, backward1x, backward1y, backward1z

    real(8) :: x_left, x_right, y_left, y_right, z_left, z_right
    real(8) :: del_x, del_y, del_z, xi, yi, zi
    real :: T1, T2, T3, T4, T5, T6, T7

    integer :: i,j,k,l, istart, iend, istart2, iend2

    if(rank==0) then
        call cpu_time(T1)
    end if

    call MPI_INIT(ierr)

    call MPI_COMM_SIZE(MPI_COMM_WORLD, nprocs, ierr)
    call MPI_COMM_RANK(MPI_COMM_WORLD, rank, ierr)

    if(rank==0) then
        print *, "Enter the domain : x_left, x_right, y_left, y_right, z_left, z_right"
        read (*,*) x_left, x_right, y_left, y_right, z_left, z_right
    end if

    call MPI_BCAST(x_left, 1, MPI_REAL8, 0, MPI_COMM_WORLD, ierr)
    call MPI_BCAST(x_right, 1, MPI_REAL8, 0, MPI_COMM_WORLD, ierr)
    call MPI_BCAST(y_left, 1, MPI_REAL8, 0, MPI_COMM_WORLD, ierr)
    call MPI_BCAST(y_right, 1, MPI_REAL8, 0, MPI_COMM_WORLD, ierr)
    call MPI_BCAST(z_left, 1, MPI_REAL8, 0, MPI_COMM_WORLD, ierr)
    call MPI_BCAST(z_right, 1, MPI_REAL8, 0, MPI_COMM_WORLD, ierr)

    call MPI_BARRIER(MPI_COMM_WORLD, ierr)

    del_x = (x_right-x_left)/nx
    del_y = (y_right-y_left)/ny
    del_z = (z_right-z_left)/nz

    ! Allocate
    if( .not. allocated(func) ) allocate(func(nx,ny,nz,3))
    if( .not. allocated(actual_curl) ) allocate(actual_curl(nx,ny,nz,3))
    if( .not. allocated(curl) ) allocate(curl(nx,ny,nz,3))
    ! if( .not. allocated(curltempx) ) allocate(curltempx(nx,ny,nz,3))
    ! if( .not. allocated(curltempy) ) allocate(curltempy(nx,ny,nz,3))
    ! if( .not. allocated(curltempz) ) allocate(curltempz(nx,ny,nz,3))
    if( .not. allocated(curltemp) ) allocate(curltemp(nx,ny,nz,9))
    
    if( .not. allocated(kappax) ) allocate(kappax(nx))
    if( .not. allocated(kappay) ) allocate(kappay(ny))
    if( .not. allocated(kappaz) ) allocate(kappaz(nz))

    if( .not. allocated(der1x) ) allocate(der1x(nx)) 
    if( .not. allocated(der1y) ) allocate(der1y(ny))
    if( .not. allocated(der1z) ) allocate(der1z(nz)) 

    if( .not. allocated(fhat1x) ) allocate(fhat1x(nx))
    if( .not. allocated(fhat1y) ) allocate(fhat1y(ny))
    if( .not. allocated(fhat1z) ) allocate(fhat1z(nz))

    if( .not. allocated(dfhat1x) ) allocate(dfhat1x(nx)) 
    if( .not. allocated(dfhat1y) ) allocate(dfhat1y(ny))
    if( .not. allocated(dfhat1z) ) allocate(dfhat1z(nz))

    ! if( .not. allocated(derx) ) allocate(derx(nx,ny,nz,3)) 
    ! if( .not. allocated(dery) ) allocate(dery(nx,ny,nz,3))
    ! if( .not. allocated(derz) ) allocate(derz(nx,ny,nz,3))
    ! derx(:,:,:,i) = der(:,:,:,i)
    ! dery(:,:,:,i) = der(:,:,:,i+3)
    ! derz(:,:,:,i) = der(:,:,:,i+6)
    if( .not. allocated(der) ) allocate(der(nx,ny,nz,9))
    der = 0
    
    if(rank == 0) then
        call cpu_time(T3)
        print *, "MPI_BCAST and variable defining", T3-T1
    end if

    call para_range(1,n,rank,istart,iend)

    do i=1,nx
        do j=1,ny
            do k=1,nz
                xi = x_left + (i-1)*del_x
                yi = y_left + (j-1)*del_y
                zi = z_left + (k-1)*del_z
                func(i,j,k,1) = (xi**2)*yi
                func(i,j,k,2) = xi*yi*zi
                func(i,j,k,3) = -1*(xi**2)*(yi**2)
                ! func(i,j,k) = xi - xi*yi + zi**2
                ! actual_curl(i,j,k,1) = 2*cos(xi)*cos(yi)*cos(zi)
                ! actual_curl(i,j,k,2) = -2*sin(xi)*cos(zi)*sin(yi)
                ! actual_curl(i,j,k,3) = -2*sin(xi)*cos(yi)*sin(zi)
            enddo
        enddo
    enddo

    if(rank==0) then
        call cpu_time(T4)
        print *, "Post assignment of func", T4-T3
    end if
    
    call dfftw_plan_dft_r2c_1d(forward1x, nx, func(:,1,1,1), fhat1x, FFTW_ESTIMATE)
    call dfftw_plan_dft_c2r_1d(backward1x, nx, dfhat1x, der1x, FFTW_ESTIMATE)
    call dfftw_plan_dft_r2c_1d(forward1y, ny, func(1,:,1,1), fhat1y, FFTW_ESTIMATE)
    call dfftw_plan_dft_c2r_1d(backward1y, ny, dfhat1y, der1y, FFTW_ESTIMATE)
    call dfftw_plan_dft_r2c_1d(forward1z, nz, func(1,1,:,1), fhat1z, FFTW_ESTIMATE)
    call dfftw_plan_dft_c2r_1d(backward1z, nz, dfhat1z, der1z, FFTW_ESTIMATE)

    call fftfreq(x_right-x_left, nx, kappax)
    call fftfreq(y_right-y_left, ny, kappay)
    call fftfreq(z_right-z_left, nz, kappaz)

    if(rank==0) then
        call cpu_time(T5)
        print *, "Calling fftw plan and finding kappa", T5-T4
    end if
    

    do i=1,n
        do j=istart, iend
            call dfftw_execute_dft_r2c(forward1x, func(:,i,j,2), fhat1x)
            dfhat1x = kappax*Z*fhat1x
            call dfftw_execute_dft_c2r(backward1x, dfhat1x, der1x)
            der(:,i,j,2) = der1x/(nx*1.0)

            call dfftw_execute_dft_r2c(forward1x, func(:,i,j,3), fhat1x)
            dfhat1x = kappax*Z*fhat1x
            call dfftw_execute_dft_c2r(backward1x, dfhat1x, der1x)
            der(:,i,j,3) = der1x/(nx*1.0)

            call dfftw_execute_dft_r2c(forward1y, func(i,:,j,1), fhat1y)
            dfhat1y = kappay*Z*fhat1y
            call dfftw_execute_dft_c2r(backward1y, dfhat1y, der1y)
            der(i,:,j,4) = der1y/(ny*1.0)

            call dfftw_execute_dft_r2c(forward1y, func(i,:,j,3), fhat1y)
            dfhat1y = kappay*Z*fhat1y
            call dfftw_execute_dft_c2r(backward1y, dfhat1y, der1y)
            der(i,:,j,6) = der1y/(ny*1.0)

            call dfftw_execute_dft_r2c(forward1z, func(i,j,:,1), fhat1z)
            dfhat1z = kappaz*Z*fhat1z
            call dfftw_execute_dft_c2r(backward1z, dfhat1z, der1z)
            der(i,j,:,7) = der1z/(nz*1.0)

            call dfftw_execute_dft_r2c(forward1z, func(i,j,:,2), fhat1z)
            dfhat1z = kappaz*Z*fhat1z
            call dfftw_execute_dft_c2r(backward1z, dfhat1z, der1z)
            der(i,j,:,8) = der1z/(nz*1.0)
        enddo
    enddo

    ! if(rank==0) then
    !     do i=1,n
    !         do j=1,n
    !             call dfftw_execute_dft_r2c(forward1x, func(:,i,j,2), fhat1x)
    !             dfhat1x = kappax*Z*fhat1x
    !             call dfftw_execute_dft_c2r(backward1x, dfhat1x, der1x)
    !             der(:,i,j,2) = der1x/(nx*1.0)
    
    !             call dfftw_execute_dft_r2c(forward1y, func(i,:,j,1), fhat1y)
    !             dfhat1y = kappay*Z*fhat1y
    !             call dfftw_execute_dft_c2r(backward1y, dfhat1y, der1y)
    !             der(i,:,j,4) = der1y/(ny*1.0)
    
    !             call dfftw_execute_dft_r2c(forward1y, func(i,:,j,3), fhat1y)
    !             dfhat1y = kappay*Z*fhat1y
    !             call dfftw_execute_dft_c2r(backward1y, dfhat1y, der1y)
    !             der(i,:,j,5) = der1y/(ny*1.0)
    !         enddo
    !     enddo
    ! end if

    ! if(rank==1) then
    !     do i=1,n
    !         do j=1,n
    !             call dfftw_execute_dft_r2c(forward1x, func(:,i,j,3), fhat1x)
    !             dfhat1x = kappax*Z*fhat1x
    !             call dfftw_execute_dft_c2r(backward1x, dfhat1x, der1x)
    !             der(:,i,j,3) = der1x/(nx*1.0)
    
    !             call dfftw_execute_dft_r2c(forward1z, func(i,j,:,1), fhat1z)
    !             dfhat1z = kappaz*Z*fhat1z
    !             call dfftw_execute_dft_c2r(backward1z, dfhat1z, der1z)
    !             der(i,j,:,7) = der1z/(nz*1.0)
    
    !             call dfftw_execute_dft_r2c(forward1z, func(i,j,:,2), fhat1z)
    !             dfhat1z = kappaz*Z*fhat1z
    !             call dfftw_execute_dft_c2r(backward1z, dfhat1z, der1z)
    !             der(i,j,:,8) = der1z/(nz*1.0)
    !         enddo
    !     enddo
    ! end if
    

    if(rank == 0) then
        call cpu_time(T6)
        print *, "Execute FFT and find derivative for one of the processors", T6-T5
    end if
    
    if(rank /= 0) then
        call MPI_SEND(istart, 1, MPI_INT, 0, rank+nprocs, MPI_COMM_WORLD, ierr)
        call MPI_SEND(iend, 1, MPI_INT, 0, rank+nprocs*2, MPI_COMM_WORLD, ierr)
        call MPI_SEND(der, n*n*n*9, MPI_REAL8, 0, rank, MPI_COMM_WORLD, ierr)
        ! call MPI_SEND(dery, n*n*n*3, MPI_REAL8, 0, rank+nprocs*3, MPI_COMM_WORLD, ierr)
        ! call MPI_SEND(derz, n*n*n*3, MPI_REAL8, 0, rank+nprocs*4, MPI_COMM_WORLD, ierr)
    end if

    if(rank==0) then
        do i=1,nprocs-1
            call MPI_RECV(istart2, 1, MPI_INT, i, i+nprocs, MPI_COMM_WORLD, status1, ierr)
            call MPI_RECV(iend2, 1, MPI_INT, i, i+nprocs*2, MPI_COMM_WORLD, status1, ierr)
            call MPI_RECV(curltemp, n*n*n*9, MPI_REAL8, i, i, MPI_COMM_WORLD, status1, ierr)
            ! call MPI_RECV(curltempy, n*n*n*3, MPI_REAL8, i, i+nprocs*3, MPI_COMM_WORLD, status1, ierr)
            ! call MPI_RECV(curltempz, n*n*n*3, MPI_REAL8, i, i+nprocs*4, MPI_COMM_WORLD, status1, ierr)
            ! Based on the received values from each of the processors, we can 
            ! rearrange them back into the master processor
            do l=1,n
                do j=istart2,iend2
                    der(:,l,j,2) = curltemp(:,l,j,2)
                    der(:,l,j,3) = curltemp(:,l,j,3)
                    der(l,:,j,4) = curltemp(l,:,j,4)
                    der(l,:,j,6) = curltemp(l,:,j,6)
                    der(l,j,:,7) = curltemp(l,j,:,7)
                    der(l,j,:,8) = curltemp(l,j,:,8)
                enddo
            enddo
        enddo

        call cpu_time(T3)
        print *, "Final send recv", T3-T6

        curl(:,:,:,1) = der(:,:,:,6) - der(:,:,:,8)
        curl(:,:,:,2) = der(:,:,:,7) - der(:,:,:,3)
        curl(:,:,:,3) = der(:,:,:,2) - der(:,:,:,4)

        ! do i=1,n
        !     do j=1,n
        !         do k=1,n
        !             print *, curl(i,j,k,:)
        !         enddo
        !     enddo
        ! enddo
        call cpu_time(T7)
        print *, "Calculating curl", T7-T3
        
    end if

    

    call dfftw_destroy_plan(forward1x, func(:,1,1,1), fhat1x)
    call dfftw_destroy_plan(backward1x, dfhat1x, der1x)
    call dfftw_destroy_plan(forward1y, func(1,:,1,1), fhat1y)
    call dfftw_destroy_plan(backward1y, dfhat1y, der1y)
    call dfftw_destroy_plan(forward1z, func(1,1,:,1), fhat1z)
    call dfftw_destroy_plan(backward1z, dfhat1z, der1z)

    if(allocated(func)) deallocate(func)
    if(allocated(actual_curl)) deallocate(actual_curl)
    if(allocated(curl)) deallocate(curl)
    ! if(allocated(curltempx)) deallocate(curltempx)
    ! if(allocated(curltempy)) deallocate(curltempy)
    ! if(allocated(curltempz)) deallocate(curltempz)
    if(allocated(curltemp)) deallocate(curltemp)

    if(allocated(der1x)) deallocate(der1x)
    if(allocated(kappax)) deallocate(kappax)
    if(allocated(der1y)) deallocate(der1y)
    if(allocated(kappay)) deallocate(kappay)
    if(allocated(der1z)) deallocate(der1z)
    if(allocated(kappaz)) deallocate(kappaz)

    if(allocated(fhat1x)) deallocate(fhat1x)
    if(allocated(fhat1y)) deallocate(fhat1y)
    if(allocated(fhat1z)) deallocate(fhat1z)

    if(allocated(dfhat1x)) deallocate(dfhat1x)
    if(allocated(dfhat1y)) deallocate(dfhat1y)
    if(allocated(dfhat1z)) deallocate(dfhat1z)

    ! if(allocated(derx)) deallocate(derx)
    ! if(allocated(dery)) deallocate(dery)
    ! if(allocated(derz)) deallocate(derz)
    if(allocated(der)) deallocate(der)

    if(rank==0) then
        call cpu_time(T2)
        print *, "Deallocate and destroy plan", T2-T7
        print *,"Completion",  T2-T1
    end if

    call MPI_FINALIZE(ierr)

end program


subroutine fftfreq(L,n,kappaf)

    implicit none
    real(8) :: L, PI = 4*atan(1.d0)
    integer :: n, i
    real(8) :: kappaf(n)
    real(8), allocatable :: kappa(:)
    ! if(.not. allocated(res)) allocate (res(n))
    if(.not. allocated(kappa)) allocate (kappa(n))

    do i=1,n
    enddo

    kappa(1) = -n/2.0
    do i=2,(n/2+1)
        kappa(i) = (i-1)-n/2.0
        kappa(n+2-i) = n/2.0-(i-1)
    enddo

    kappaf(1:n/2) = kappa(CEILING(n/2.0)+1:n)*2*PI/L
    kappaf(n/2+1:n) = kappa(1:CEILING(n/2.0)+1)*2*PI/L

    ! print *, res

    if(allocated(kappa)) deallocate(kappa)
    ! if(allocated(res)) deallocate(res)
    return

end subroutine fftfreq


subroutine para_range(n1, n2, irank, istart, iend)
    
    use nodeinfo

    implicit none
    integer :: n1, n2, irank, istart, iend, iwork1, iwork2
    
    iwork1 = (n2-n1+1)/nprocs
    iwork2 = MOD(n2-n1+1, nprocs)

    istart = n1 + irank*iwork1 + MIN(irank, iwork2)
    iend = istart + iwork1 - 1

    if(iwork2 > irank) iend = iend+1
    return

end subroutine para_range
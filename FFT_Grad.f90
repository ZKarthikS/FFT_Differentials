! Program to find the Gradient of a 3 dimensional scalar function by serially computing the FFT

program fft_grad

    implicit none

    include '/usr/local/include/fftw3.f'

    integer, parameter :: n = 250
    integer, parameter :: nx = n, ny = n, nz = n
    complex, parameter :: Z = (0,1)

    double complex, allocatable :: fhat1x(:), fhat1y(:), fhat1z(:), dfhat1x(:), dfhat1y(:), dfhat1z(:)

    real(8), allocatable :: actual_grad(:,:,:,:), grad(:,:,:,:), func(:,:,:)
    real(8), allocatable :: kappax(:), kappay(:), kappaz(:)
    real(8), allocatable :: der1x(:), der1y(:), der1z(:)

    integer(8) :: forward1x, forward1y, forward1z, backward1x, backward1y, backward1z

    real(8) :: x_left, x_right, y_left, y_right, z_left, z_right
    real(8) :: del_x, del_y, del_z, xi, yi, zi
    real :: T1, T2, T3, T4, T5

    integer :: i,j,k

    call cpu_time(T1)

    ! print *, "Enter the domain : x_left, x_right, y_left, y_right, z_left, z_right"
    read (*,*) x_left, x_right, y_left, y_right, z_left, z_right

    del_x = (x_right-x_left)/nx
    del_y = (y_right-y_left)/ny
    del_z = (z_right-z_left)/nz

    ! Allocate
    if( .not. allocated(func) ) allocate(func(nx,ny,nz))
    if( .not. allocated(actual_grad) ) allocate(actual_grad(nx,ny,nz,3))
    if( .not. allocated(grad) ) allocate(grad(nx,ny,nz,3))
    
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
    
    call cpu_time(T3)
    ! print *, "Before beginning", T3-T1

    do i=1,nx
        do j=1,ny
            do k=1,nz
                xi = x_left + (i-1)*del_x
                yi = y_left + (j-1)*del_y
                zi = z_left + (k-1)*del_z
                func(i,j,k) = 2*sin(xi)*cos(yi)*cos(zi)
                ! func(i,j,k) = xi - xi*yi + zi**2
                ! actual_grad(i,j,k,1) = 2*cos(xi)*cos(yi)*cos(zi)
                ! actual_grad(i,j,k,2) = -2*sin(xi)*cos(zi)*sin(yi)
                ! actual_grad(i,j,k,3) = -2*sin(xi)*cos(yi)*sin(zi)
            enddo
        enddo
    enddo

    ! do i=1,nx
    !     do j=1,ny
    !         print *, func(i,j,:)
    !     enddo
    ! enddo

    call cpu_time(T4)
    ! print *, "Post assignment of func", T4-T3

    call dfftw_plan_dft_r2c_1d(forward1x, nx, func(:,1,1), fhat1x, FFTW_ESTIMATE)
    call dfftw_plan_dft_c2r_1d(backward1x, nx, dfhat1x, der1x, FFTW_ESTIMATE)
    call dfftw_plan_dft_r2c_1d(forward1y, ny, func(1,:,1), fhat1y, FFTW_ESTIMATE)
    call dfftw_plan_dft_c2r_1d(backward1y, ny, dfhat1y, der1y, FFTW_ESTIMATE)
    call dfftw_plan_dft_r2c_1d(forward1z, nz, func(1,1,:), fhat1z, FFTW_ESTIMATE)
    call dfftw_plan_dft_c2r_1d(backward1z, nz, dfhat1z, der1z, FFTW_ESTIMATE)

    call fftfreq(x_right-x_left, nx, kappax)
    call fftfreq(y_right-y_left, ny, kappay)
    call fftfreq(z_right-z_left, nz, kappaz)

    ! do j=1,ny
    !     do k=1,nz            
    !         call dfftw_execute_dft_r2c(forward1x, func(:,j,k), fhat1x)

    !         dfhat1x = kappax*Z*fhat1x

    !         call dfftw_execute_dft_c2r(backward1x, dfhat1x, der1x)

    !         grad(:,j,k,1) = der1x/(nx*1.0)
    !     enddo
    ! enddo

    ! do i=1,nx
    !     do k=1,nz
    !         call dfftw_execute_dft_r2c(forward1y, func(i,:,k), fhat1y)

    !         dfhat1y = kappay*Z*fhat1y

    !         call dfftw_execute_dft_c2r(backward1y, dfhat1y, der1y)

    !         grad(i,:,k,2) = der1y/(ny*1.0)
    !     enddo
    !     do j=1,ny
    !         call dfftw_execute_dft_r2c(forward1z, func(i,j,:), fhat1z)

    !         dfhat1z = kappaz*Z*fhat1z

    !         call dfftw_execute_dft_c2r(backward1z, dfhat1z, der1z)

    !         grad(i,j,:,3) = der1z/(nz*1.0)
    !     enddo
    ! enddo

    call cpu_time(T5)
    ! print *, "Calling fftw plan and finding kappa", T5-T4

    do i=1,n
        do j=1,n
            call dfftw_execute_dft_r2c(forward1x, func(:,i,j), fhat1x)
            dfhat1x = kappax*Z*fhat1x
            call dfftw_execute_dft_c2r(backward1x, dfhat1x, der1x)
            grad(:,i,j,1) = der1x/(nx*1.0)

            call dfftw_execute_dft_r2c(forward1y, func(i,:,j), fhat1y)
            dfhat1y = kappay*Z*fhat1y
            call dfftw_execute_dft_c2r(backward1y, dfhat1y, der1y)
            grad(i,:,j,2) = der1y/(ny*1.0)

            call dfftw_execute_dft_r2c(forward1z, func(i,j,:), fhat1z)
            dfhat1z = kappaz*Z*fhat1z
            call dfftw_execute_dft_c2r(backward1z, dfhat1z, der1z)
            grad(i,j,:,3) = der1z/(nz*1.0)
        enddo
    enddo

    do i=1,nx
        do j=1,ny
            print *, grad(i,j,:,3)
        enddo
    enddo

    

    call dfftw_destroy_plan(forward1x, func(:,1,1), fhat1x)
    call dfftw_destroy_plan(backward1x, dfhat1x, der1x)
    call dfftw_destroy_plan(forward1y, func(1,:,1), fhat1y)
    call dfftw_destroy_plan(backward1y, dfhat1y, der1y)
    call dfftw_destroy_plan(forward1z, func(1,1,:), fhat1z)
    call dfftw_destroy_plan(backward1z, dfhat1z, der1z)

    if(allocated(func)) deallocate(func)
    if(allocated(actual_grad)) deallocate(actual_grad)
    if(allocated(grad)) deallocate(grad)

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

    call cpu_time(T2)

    ! print *, "Calculating grad and destroy plan", T2-T5

    ! print *,"Completion",  T2-T1

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
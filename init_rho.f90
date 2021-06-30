subroutine init_rho(rho, n, rho_low, rho_high)
    implicit none
    real, dimension(n) ::rho
    real rho_low, rho_high
    real dr
    integer i,n

!   do i=1,2*n/6
!       rho(i)=rho_high
!   end do
!   dr = log10(rho_low/rho_high) *3.0/n
!   do i=1+2*n/6, 3*n/6
!       rho(i) = rho_high * 10**((i-2*n/6)*dr)
!   end do
!   do i=1+3*n/6,n
!       rho(i)=rho_low
!   end do
    
    do i=1, n/2
        rho(i) = rho_high
    end do
    do i=1+n/2, n
        rho(i) = rho_low
    end do


    return

    end

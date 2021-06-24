subroutine init_rho(rho, n, rho_low, rho_high)
    implicit none
    real, dimension(n) ::rho
    real rho_low, rho_high
    integer i,n

    do i=1,n/2
        rho(i)=rho_high
    end do
    do i=1+n/2,n
        rho(i)=rho_low
    end do


    return

    end

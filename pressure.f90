real function Press(rho, A, gam)
    implicit none
    real :: rho, A, gam

    Press = A * rho**gam

    end function Press

subroutine Pressure_array(A, gam, rho, P, n)
    implicit none
    integer :: i,n
    real :: A, gam
    real, dimension(n) :: rho, P
    real :: Press

    do i=1,n
        P(i) = Press(rho(i), A, gam)
    end do

    return
    end

subroutine pres_diff(P, dP, n, dx, rho_low, rho_high, A, gam)
    implicit none
    real, dimension(n) :: P, dP
    real :: Press
    real :: dx, dx2, rho_low, rho_high, gam, A
    integer :: i,n

    dx2 = 2*dx
    do i=2, n-1
        dP(i) = (P(i-1) - P(i+1))/dx2
    end do
    dP(1) = (Press(rho_low, A, gam) - P(2))/dx2
    dP(n) = (P(n-1) - Press(rho_high, A, gam))/dx2
    return
    end 

subroutine Sound_Speed(P, rho, c, n, gam)
    IMPLICIT none
    integer :: i, n
    real, dimension(n) :: P, rho, c
    real :: gam

    do i=1,n
        c = sqrt(gam * P(i)/rho(i))
    end do

    RETURN
    end

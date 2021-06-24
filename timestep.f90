SUBROUTINE Get_timestep(c,v,n, eps, dx, dt)
    IMPLICIT none
    integer :: n,i
    REAL, dimension(n) :: c,v
    REAL :: eps, dx
    REAL :: dt, dt_1, maxv

    dt_1 = 1.0/max(c(1), v(1))
    do i=2,n
        maxv = 1.0/max(c(i), v(i))
        dt_1 = min(maxv, dt_1)
    END do

    dt = eps*dx/dt_1

    RETURN
    end
   

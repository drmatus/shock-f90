SUBROUTINE timestep(c,v,n, eps, dx, dt)
    IMPLICIT none
    integer :: n,i
    REAL, dimension(n) :: c,v
    REAL :: eps, dx
    REAL :: dt, dt_1, maxv

    dt_1 = 1.0/max(abs(c(1)), abs(v(1)))
!   dt_1 = 1.0/(abs(c(1))+ abs(v(1)))
    do i=2,n
!       maxv = 1.0/(abs(c(i))+ abs(v(i)))
        maxv = 1.0/max(abs(c(i)), abs(v(i)))
        dt_1 = min(maxv, dt_1)
    END do

    dt = eps*dx*dt_1
    if (dt < 100.0) then
        dt = 100.0
    end if

    RETURN
    end
   

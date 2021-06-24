subroutine init_x(x, n, L, dx)
    implicit none
    real, dimension(n) ::x
    real L, dx
    integer i,n

    dx = L/n

    do i=1,n
        x(i) = dx*i
    end do

    return

    end

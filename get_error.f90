subroutine get_error(x1, x2, n, error)
    implicit none
    real, dimension(n) :: x1, x2
    real :: error, tmp
    integer :: n, i

    error = 0.0
    do i =1,n
        tmp = abs(x1(i) - x2(i))
        if (x1(i).ne.0.0) then
            tmp = tmp/abs(x1(i))
        end if 
        error = max(tmp, error)
    end do


    return
    end

subroutine mk_matrix_B(B, v, n,  alpha, theta)
    implicit none
    real, dimension(n,n) :: B
    real, dimension(n)   :: v
    real                 :: alpha, theta
    integer              :: n, i

    B(1,1) = 1
    B(1,2) = -alpha*v(2)
    B(n,n) = 1
    B(n,n-1) = alpha*v(n-1)

    do i=2, n-1
        B(i, i-1) = (1-theta) * alpha * v(i-1)
        B(i,i)    =  1
        B(i, i+1) = (theta-1) * alpha * v(i+1)
    end do

    return
    end

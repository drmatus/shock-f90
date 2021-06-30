subroutine mk_matrix_A(Aa, Ab, Ac, v1, n,  alpha, theta)
    implicit none
    real, dimension(n) :: Aa, Ab, Ac
    real, dimension(n)   :: v1
    real                 :: alpha, theta
    integer              :: n, i

    Ab(1) = 1.0
    Ac(1) = theta*alpha*v1(2)
    Aa(n) = -theta*alpha*v1(n-1)
    Ab(n) = 1.0
    do i=2,n-1
        Aa(i) = -theta*alpha*v1(i-1)
        Ab(i) = 1.0
        Ac(i) = theta*alpha*v1(i+1)
    end do

    return
    end

subroutine mk_matrix_A(A, v1, n,  alpha, theta)
    implicit none
    real, dimension(n,n) :: A
    real, dimension(n)   :: v1
    real                 :: alpha, theta
    integer              :: n, i

    A(1,1)   = 1
    A(1,2)   =  theta*alpha*v1(2)
    A(n,n)   = 1
    A(n,n-1) = -theta*alpha*v1(n-1)

    do i=2, n-1
        A(i, i-1) = -theta * alpha * v1(i-1)
        A(i,i)    =  1
        A(i, i+1) =  theta * alpha * v1(i+1)
    end do

    return
    end

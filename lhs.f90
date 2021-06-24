SUBROUTINE lhs(A,x,y,n)
    IMPLICIT none
    REAL, dimension(n,n) :: A
    REAL, dimension(n)   :: x,y
    INTEGER              :: i,n

    y = 0.0
    DO i=1,n
        y = y + x(i) * A(:,i) 
    END DO
    RETURN
    END 

subroutine save_timestep(x, rho, v, P, n, t, nstep)
    implicit none
    real, dimension(n) :: x, rho, v, P
    real               :: t
    integer            :: n, nstep, i
    character*80       :: nombre

!   Create filename
    write(nombre, 100) "out_",nstep,".dat"
!   Open new file
    open(10, FILE=nombre, ACTION='WRITE')

    write(10, 200) n, t, nstep
    do i=1,n
        write(10, 300) x(i), rho(i), v(i), P(i)
    end do

100 FORMAT(A4,i5.5,A4)
200 FORMAT(I3, F10.3, I3)
300 FORMAT(4(ES14.7E2))

    return 
    end

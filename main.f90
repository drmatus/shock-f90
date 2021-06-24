PROGRAM main
    implicit none
    real,dimension(:), allocatable ::  rho, v, phi
    real,dimension(:), allocatable ::  rho1, v1, phi1
    real,dimension(:), allocatable ::  x, c
    real,dimension(:), allocatable ::  P, dP, P1, dP1
    REAL, dimension(:,:), allocatable :: A_coef, B_coef
    REAL, dimension(3,3) :: m_test
    REAL, dimension(3) :: x_test, y_test

    real :: rho_low, rho_high
    real :: dx, dt, alpha
    REAL :: t, t_end
    real :: L, A, gam
    REAL :: eps, theta
    integer :: i, n, step

    PRINT *, "== init test =="
    DO i=1,3
        x_test(i) = i
        DO step=1,3
            m_test(step,i) = i*step
        END DO
    END DO

    PRINT *, x_test
    PRINT *, m_test
    CALL lhs(m_test, x_test, y_test, 3)
    PRINT *, y_test
    PRINT *, "== Fin test =="

!   Set the number of points of the grid
    n = 30
    rho_high = 1e-20
    rho_low = rho_high*1e-4
    L = 1e15
    gam = 5.0/3.0
    A = 1.78e22 ! Del informe
    t = 0
    t_end = 1e17
    eps = 0.1
    theta = 0.5 !0 = explicito, 1= implicito
    step = 0

!   Allocate memory to the arrays
    allocate(x(n))
    allocate(c(n))
    allocate(rho(n))
    allocate(rho1(n))
    allocate(v(n))
    allocate(v1(n))
    allocate(phi(n))
    allocate(phi1(n))
    allocate(P(n))
    allocate(dP(n))
    allocate(P1(n))
    allocate(dP1(n))
    allocate(A_coef(n,n))
    allocate(B_coef(n,n))

!   Initialize arrays with some values
    v = 0.0
    v1 = 0.0
    rho1 = 0.0
    phi = 0.0
    phi1 = 0.0
    A_coef = 0.0
    B_coef = 0.0


    call init_x(x,n,L,dx)
    call init_rho(rho,n, rho_low, rho_high)
    print *,c

    print *, "start sim loop"
    DO while (t < t_end)

        CALL Sound_Speed(P, rho, c, n, gam)
        CALL Get_timestep(c,v,n,eps,dx,dt)
        print *, step, dt

        alpha = dt/(2*dx)

        !======= Paso Predictor =======
        ! El paso predictor consiste en el metodo de eluler explicito.
        ! Este paso solo requiere de la matriz B, por lo que construimos dicha 
        ! matriz usando theta = 0.
        CALL mk_matrix_B(B_coef, v , n, alpha, 0)

        ! Para el calculo del momento, se necesita la presion, que es funcion de rho, y su derivada.
        CALL pressure_array(A, gam, rho, P, n)
        CALL pres_diff(P, dP, n, dx, rho_low, rho_high, A, gam)

        !Para la conservacion de masa, multiplicamos la matriz B por rho
        Call lhs(B_coef, rho, rho1, n)

        ! Para la conservacion de momento, hacemos lo mismo, y sumamos la derivada de la presion:
        Call lhs(B_coef, phi, phi1, n)



        !======= Fin Paso Predictor =======

        call save_timestep(x, rho, v,P,n, t, step)

    t = t + dt
    step = step +1

    END DO
    print *, "END sim loop"

    call pres_diff(P, dP, n, dx, rho_low, rho_high, A, gam)


! ======== FIN DEL PROGRAMA ===========
    deallocate(x)
    deallocate(rho)
    deallocate(rho1)
    deallocate(v)
    deallocate(v1)

END PROGRAM

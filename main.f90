PROGRAM main
    implicit none
    real,dimension(:), allocatable ::  rho, v, phi
    real,dimension(:), allocatable ::  rho1, v1, phi1
    real,dimension(:), allocatable ::  x, c
    real,dimension(:), allocatable ::  P, dP, P1, dP1
    REAL, dimension(:,:), allocatable :: A_coef, B_coef

    real :: rho_low, rho_high
    real :: dx, dt, alpha
    REAL :: t, t_end
    real :: L, A, gam
    REAL :: eps, theta
    integer :: i, n, step

!   Set the number of points of the grid
    n = 30
    rho_high = 1e-20
    rho_low = rho_high*1e-4
    L = 1e15
    gam = 5.0/3.0
    A = 1.78e2 ! Del informe
    t = 0
    t_end = 1e7
    eps = 0.01
    theta = 0.0 !0 = explicito, 1= implicito
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

    print *, "start sim loop"
    DO while (t < t_end)

        CALL Sound_Speed(P, rho, c, n, gam)
        CALL timestep(c,v,n,eps,dx,dt)
        print *,c
        print *, step, dt

        alpha = dt/(2*dx)

        !======= Paso Predictor =======
        ! El paso predictor consiste en el metodo de euler explicito.
        ! Este paso solo requiere de la matriz B, por lo que construimos dicha 
        ! matriz usando theta = 0.
        CALL mk_matrix_B(B_coef, v , n, alpha, 0)

        ! Para el calculo del momento, se necesita la presion, que es funcion de rho, y su derivada.
        CALL pressure_array(A, gam, rho, P, n)
        CALL pres_diff(P, dP, n, dx, rho_low, rho_high, A, gam)

        !Para la conservacion de masa, multiplicamos la matriz B por rho
        Call lhs(B_coef, rho, rho1, n)

        ! Para la conservacion de momento, hacemos lo mismo, y sumamos la derivada de la presion:
        do i=1,n
            phi(i) = rho(i)*v(i)
        end do
        
        Call lhs(B_coef, phi, phi1, n)
        phi1 = phi1 - dt*dP 

        ! Obrenemos v a partir de phi y rho:
        DO i=1,n
            v1(i) = phi1(i)/rho1(i)
        END DO
        !======= Fin Paso Predictor =======

        ! Si solo queremos el metodo de euler explicito, entonces solo necesitamos salta al final del timestep:
        IF (theta.EQ.0) THEN
            goto 100
        END IF 

        !===== Pasos correctores =====

        

        !===== Fin Pasos correctores =====

        ! Una vez conformes con los valores rho1 y v1, reemplazamos rho y v:
100     DO i=1,n
            rho(i) = rho1(i)
            v(i)   = v1(i)
        END DO
        ! Y guardamos los resultados.
        call save_timestep(x, rho, v,P, n, t, step)

    t = t + dt
    step = step +1

    END DO
    print *, "END sim loop"

! ======== FIN DEL PROGRAMA ===========
    deallocate(x)
    deallocate(rho)
    deallocate(rho1)
    deallocate(v)
    deallocate(v1)

END PROGRAM

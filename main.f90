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

!   ======= Program parameters ==========
!   Set the number of points of the grid
    n = 10
!   Densitty of the first half
    rho_high = 1e-20
!   Densitty of the second half
    rho_low = rho_high*1e-4
!   Lenght of simulation area
    L = 3e15
!   Adiabatic constant
    gam = 5.0/3.0
!   time at end of simulation
    t_end = 1e30
!   staility parameter
    eps = 0.00001
!   scheme selector
    theta = 0.0 !0 = explicito, 1= implicito


!   ======= Other variable initialization ==========
    alpha = 0.0
    dt = 0.0
    dx = 0.0
    step = 0.0
    t = 0.0
    i = 0


!   ======= Allocate memory to the arrays ========
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

!   ========== Initialize arrays with some values ========
    v = 0.0
    v1 = 0.0
    rho1 = 0.0
    phi = 0.0
    phi1 = 0.0
    A_coef = 0.0
    B_coef = 0.0

!   Set the normalization constant for the gas:
    CALL gas_norm(rho_high, gam, 10.0, A)

    call init_x(x,n,L,dx)
    call init_rho(rho,n, rho_low, rho_high)
    CALL pressure_array(A, gam, rho, P, n)

    call save_timestep(x, rho, v, P, n, t, step)
    print *, step, t

    print *, "start sim loop"
    DO while (t < t_end)

        CALL Sound_Speed(rho, c, n, gam,A)
        CALL timestep(c,v,n,eps,dx,dt)

        alpha = dt/(2*dx)

        !======= Paso Predictor =======
        ! El paso predictor consiste en el metodo de euler explicito.
        ! Este paso solo requiere de la matriz B, por lo que construimos dicha 
        ! matriz usando theta = 0.
        CALL mk_matrix_B(B_coef, v , n, alpha, 0)

        ! Para el calculo del momento, se necesita la presion, que es funcion de rho, y su derivada.
        CALL pressure_array(A, gam, rho, P, n)
        CALL pres_diff(P, dP, n, dx, rho_low, rho_high, A, gam)

        ! Consruimos el momentum del timestep actual
        do i=1,n
            phi(i) = rho(i)*v(i)
        end do

        !Para la conservacion de masa, multiplicamos la matriz B por rho
        Call lhs(B_coef, rho, rho1, n)

        ! Para la conservacion de momento, hacemos lo mismo, y sumamos la derivada de la presion:
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
        ! Reclculamos la presion
        CALL pressure_array(A, gam, rho, P, n)

        ! Y guardamos los resultados.

        t = t + dt
        step = step +1

        print *, step, t, dt
        call save_timestep(x, rho, v,P, n, t, step)

    END DO
    print *, "END sim loop"

! ======== FIN DEL PROGRAMA ===========
    deallocate(x)
    deallocate(rho)
    deallocate(rho1)
    deallocate(v)
    deallocate(v1)

END PROGRAM

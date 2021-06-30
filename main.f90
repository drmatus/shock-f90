PROGRAM main
    implicit none
    real,dimension(:), allocatable ::  rho, v, phi
    real,dimension(:), allocatable ::  rho1, v1, phi1
    real,dimension(:), allocatable ::  rho2, v2, phi2
    real,dimension(:), allocatable ::  x, c
    real,dimension(:), allocatable ::  P, dP, P1, dP1
    REAL, dimension(:,:), allocatable :: B_coef
    real, dimension(:), allocatable :: A_coefa, A_coefb, A_coefc

    real :: rho_low, rho_high
    real :: dx, dt, alpha
    REAL :: t, t_end
    real :: L, A, gam
    REAL :: eps, theta
    real :: verr, dv
    real :: rerr, dr
    integer :: i, n, step, psteps

!   ======= Program parameters ==========
!   Set the number of points of the grid
    n = 30
!   Densitty of the first half
    rho_high = 1e-16
!   Densitty of the second half
    rho_low = rho_high*1e-4
!   Lenght of simulation area
    L = 3e15
!   Adiabatic constant
    gam = 5.0/3.0
!   time at end of simulation
    t_end = 1e11
!   staility parameter
    eps = 1.0e-1
!   scheme selector
    theta = 0.0 !0 = explicito, 1= implicito
!   relative error for rho and v:
    verr = 0.01
    rerr = 0.01


!   ======= Other variable initialization ==========
    alpha = 0.0
    dt = 0.0
    dx = 0.0
    step = 0
    psteps = 0
    t = 0.0
    i = 0
    

!   ======= Allocate memory to the arrays ========
    allocate(x(n))
    allocate(c(n))
    allocate(rho(n))
    allocate(rho1(n))
    allocate(rho2(n))
    allocate(v(n))
    allocate(v1(n))
    allocate(v2(n))
    allocate(phi(n))
    allocate(phi1(n))
    allocate(phi2(n))
    allocate(P(n))
    allocate(dP(n))
    allocate(P1(n))
    allocate(dP1(n))
    allocate(A_coefa(n))
    allocate(A_coefb(n))
    allocate(A_coefc(n))
    allocate(B_coef(n,n))

!   ========== Initialize arrays with some values ========
    v = 0.0
    v1 = 0.0
    v2 = 0.0
    rho1 = 0.0
    rho2 = 0.0
    phi = 0.0
    phi1 = 0.0
    phi2 = 0.0
    A_coefa = 0.0
    A_coefb = 0.0
    A_coefc = 0.0
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

        dv = 10.0
        dr = 10.0
        psteps = 0
        CALL Sound_Speed(rho, c, n, gam,A)
        CALL timestep(c,v,n,eps,dx,dt)
        ! dt =  3711125.750

        alpha = dt/(2*dx)

        !======= Paso Predictor =======
        ! El paso predictor consiste en el metodo de euler explicito.
        ! Este paso solo requiere de la matriz B, por lo que construimos dicha 
        ! matriz usando theta = 0.
        CALL mk_matrix_B(B_coef, v , n, alpha, 0.0)

        ! Para el calculo del momento, se necesita la presion, que es funcion de rho, y su derivada.
        CALL pressure_array(A, gam, rho, P, n)
        CALL pres_diff(P, dP, n, dx, rho_low, rho_high, A, gam)

        ! Consruimos el momentum del timestep actual
        do i=1,n
            phi(i) = rho(i)*v(i)
        end do

        !Para la conservacion de masa, multiplicamos la matriz B por rho, guardamos en rho1
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

        do while ((rerr.le.dr).or.(verr.le.dv))

        !Calcular matriz de coeficientes A, usando v1:
            call mk_matrix_A(A_coefa, A_coefb, A_coefc, v1, n, alpha, theta)

        !CalcAlar matriz de coeficientes B, usando v:
            call mk_matrix_B(B_coef, v , n, alpha, theta)

        !Calcular P1 y dP1, a partir de rho1:
            call pressure_array(A, gam, rho1, P1, n)
            call pres_diff(P1, dP1, n, dx, rho_low, rho_high, A, gam)

        !Resolver los sistemas de ecuaciones:

            !Para rho:
            !    A_coef rho2 = B_coef rho
            call lhs(B_coef, rho, rho2, n)
            ! tridag salio del numerical recipes
            call tridag(A_coefa, A_coefb, A_coefc, rho2, rho2, n)

            !Para phi:
            !    A_coef phi2 = B_coef phi - dt((1-theta)dP + theta*dP1)
            call lhs(B_coef, phi, phi2, n)
            phi2 = phi2 - dt*((1-theta)*dP  + theta*dP1)
            call tridag(A_coefa, A_coefb, A_coefc, phi2, phi2, n)
            
            call get_error(rho1, rho2, n, dr)
            call get_error(phi1, phi2, n, dv)
            rho1 = rho2
            phi1 = phi2
            DO i=1,n
                v1(i) = phi1(i)/rho1(i)
            END DO
            psteps = psteps + 1
            print *, "predictor step: ", psteps
        end do

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

        if (step > 1000) then
            exit
        end if

    END DO
    print *, "END sim loop"

! ======== FIN DEL PROGRAMA ===========

    deallocate(x)
    deallocate(c)
    deallocate(rho)
    deallocate(rho1)
    deallocate(rho2)
    deallocate(v)
    deallocate(v1)
    deallocate(v2)
    deallocate(phi)
    deallocate(phi1)
    deallocate(phi2)
    deallocate(P)
    deallocate(dP)
    deallocate(P1)
    deallocate(dP1)
    deallocate(A_coefa)
    deallocate(A_coefb)
    deallocate(A_coefc)
    deallocate(B_coef)

END PROGRAM

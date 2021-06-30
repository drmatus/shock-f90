subroutine gas_norm(rho, gam, gas_temp, A)
    implicit none
    real :: rho, gas_temp
    real :: gam, k_b, m
    real :: A

    !boltzmann constant:
    k_b =  1.38e-16  !erg K^-1
    ! Masa Hydrogeno:
    m   = 1.67e-24   ! gr

    !Usando rho en g cm^-3, temperatura en kelvin, calcular la constante A:
    ! Usando la ecuacion de estado para un gas adiavatico:
    !      P = A*rho^gamma
    ! Y la ecuacion de estado para un gas ideal:
    !      P = rho k_b*T/m
    ! Se igualan los valores P y despeja la constante A:
    A = rho**(1-gam) * gas_temp * k_b/m


    return 
    end subroutine

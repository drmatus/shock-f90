program Proyecto_B
	
	implicit none
	integer, parameter :: NX = 10, NT = 100								! Number of data points in space grid (NX) and time grid (NT).
	double precision :: XF = 5.d+0, TF = 40.d+0							! Right boundary condition. Value for domain right limit.
	double precision, dimension (NX+1) :: XD							! Grid for space data points.
	double precision, dimension (NT+1) :: TD							! Grid for time data points.

	call Space(XF, NX, XD)
	call Tiempo(TF, NT, TD)
	
	print *, XD
	print "(/)"
	print *, TD
	
	contains
		
		! This subroutine assings values for the *num* data points between 0 and
		! *rightlim*, and returns the grid (spatial discretization).
		
		subroutine Space(rightlim, num, grid)							
			implicit none
			integer :: i												! Counter to add values to the grid.																		
			integer, intent(in) :: num									! Number of data points.								
			double precision, intent(in) :: rightlim					! Right boundary condition.
			double precision :: delta									! Step.
			double precision, dimension(num+1), intent(out) :: grid		! Grid array.
			delta = rightlim/num										! Calculates the step		
			grid(1) = 0.d+0												! Left boundary condition. Sets minimum value of the domain.
			grid(num+1) = rightlim										! Right boundary condition. Sets maximum value of the domain.
			do i = 2, num												
				grid(i) = (i-1) * delta									! Fills the values in between the two extremes.
			end do
		end subroutine Space
		
		! This subroutine assings values for the *num* data points between 0 and
		! *rightlim*, and returns the grid (Time Discretization).
		
		subroutine Tiempo(tmax, num, time)								! Same behaviour as above but for time grid.
			implicit none
			integer :: i
			integer, intent(in) :: num
			double precision, intent(in) :: tmax
			double precision :: delta
			double precision, dimension(num+1), intent(out) :: time
			delta = tmax/num
			time(1) = 0.d+0
			time(num+1) = tmax
			do i = 2, num
				time(i) = (i-1) * delta
			end do
		end subroutine Tiempo

end program Proyecto_B

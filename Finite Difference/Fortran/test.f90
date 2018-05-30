program test

  implicit none

    ! Declaring Variables
    integer L, N, tf, i
    real dx, dt, x, t0


    ! Quantifying
    L = 20
    N = 201
    tf = 1
    dx = L/(N-1)
    dt = .004
	t0 = 0.1
	x = 20


    call concentration(x,t0)

	!do j = 1
    !end do

end program test

! Function for analytical concentration
subroutine concentration(x,t)

    ! Declares variables
    real x0,t0,c
	
    ! Quantifying initial conditions
    x0 = 10;
    t0 = 0.1;

    ! Analytical Function
    c = sqrt(t0/t)* (-(x-x0)**2/(4*t))

    ! Spits out concentration
    print*, c
	
	print*, 'hello'
	
	return


end subroutine concentration

	

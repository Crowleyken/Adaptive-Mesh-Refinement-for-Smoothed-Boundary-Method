! A fortran95 program for G95

program main

  implicit none

    ! Declaring Variables
    integer L, N, tf
    real dx, dt

    ! Quantifying
    L = 20
    N = 201
    tf = 1
    dx = L/(N-1)
    dt = .004

    ! Declares x vector
    real, dimension(N) :: x

    ! Declares initial x array
    x(1) = 0

    ! Loop that creates equal spacing for x
    do i = 1,N

        ! adds dt every time
        x(i+1) = x(i) + dt

    end do

    call concentration(x,t0)

   ! do j = 1



    !end do

end main


! Function for analytical concentration
subroutine concentration(x,t)

    ! Declares variables
    real c,x0,t0

    ! Quantifying initial conditions
    x0 = 10;
    t0 = 0.1;

    ! Analytical Function
    c = sqrt(t0/t)*exp(-(x-x0)^2/(4*t))

    ! Spits out concentration
    return c

end subroutine concentration

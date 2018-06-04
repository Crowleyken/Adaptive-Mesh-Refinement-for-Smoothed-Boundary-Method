! This program solves the 1D problem from Hw3 pdf and writes the 2D array to a data file

program main
  
  use Arrays
  implicit none
  
	! Interface to allow subroutine to be used
	interface
		subroutine concentration(x,t,c,x0,t0)
		! Declares variables in parameter
		real*8, intent(in), dimension(:) :: x
		real*8, intent(in) :: t, x0, t0
		real*8, intent(out), dimension(:) :: c
		end subroutine concentration
	end interface
  
  
    ! Declaring Variables
    integer L, N, tf, time, k, j
    real*8 dx, dt, t0, x0
	
	! Allowing matrices to allocate data later on, x and c are 1D, cstore is 2D
	real*8, allocatable :: x(:), c(:), cstore(:,:)
	
    ! Quantifying
    L = 20
    N = 201
	
	t0 = 0.1
    tf = 1
	
	x0 = 10
    dx = .1
    dt = .004
	
	! amount of time steps, floor rounds up to nearest integer for do loop to work down below
	time = floor((tf-t0)/dt)

	! Allocates 1D array size for x and c (from 1 to N)
	allocate(x(1:N))
	allocate(c(1:N))
	
	! Allocated 2D array size for each concentratrion gradient to be stored, from 1 to N and 1 to time
	allocate(cstore(1:N,1:time))

    ! Loop that creates equal spacing (by dx) from 1 to N. Creates x array
    do k = 1,N
        ! adds dx every time
        x(k) = dx*(k-1)
    end do
	
    
	! Initial c is found before the do loop
	call concentration(x,t0,c,x0,t0)
	
	
	! Loops from 1 to N, determining new concentration gradient at each time step
	do j = 1,N
	
	! Finite difference equation
	c = c + (dt/dx**2)*(cshift(c,-1)-2*c+cshift(c,1))
	
	! Stores each concentration gradient from each time step into a 2D array
	cstore(j,:) = c
	
	end do
	
	! calls the array module, writing a two dimensional array to a data file
	call TwoD_write(cstore,N,time)
	
end program main


! Function for analytical concentration
subroutine concentration(x,t,c,x0,t0)

	implicit none
	
    ! Declares variables in parameter
	real*8, intent(in), dimension(:) :: x
	real*8, intent(in) :: t, x0, t0
	real*8, intent(out), dimension(:) :: c
	
    ! Analytical Function
    c = sqrt(t0/t)*exp(-(x-x0)**2/(4*t))
end subroutine concentration




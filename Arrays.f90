! This module writes different dimension array sizes to a data file.

module Arrays
contains
	! Writes a 1D array to a data file
	subroutine OneD_write(Array, N1)
	implicit none
	
	real*8, intent(in), dimension(:) :: Array
	integer, intent(in) :: N1
	
	
	! opens unit, taking in the FileName parameter, unformatted (data file), replaces if there is already a data file, and writes to a data file
	open(unit=1, file='/Users/Kendell/Desktop/School/MSU/concdata.dat', form='unformatted', status='replace', action='write')
	
	! Writes the length of the columns in the 1D array (row vector)
	write(1, *) N1
	
	! Writes entire array up to its length
	write(1, *) Array(1:N1)
	
	close(1)
	end subroutine OneD_write
	
	
	
	! Writes a 2D array to a data file
	subroutine TwoD_write(Array, N1, N2)
	implicit none
	
	real*8, intent(in), dimension(:,:) :: Array
	integer, intent(in) :: N1, N2
	
	! opens unit, taking in the FileName parameter, unformatted (data file), replaces if there is already a data file, and writes to a data file
	open(unit=1, file='/Users/Kendell/Desktop/School/MSU/concdata.dat', form='unformatted', status='replace', action='write')
	
	! Writes the length of the N1 in the 2D array
	write(1, *) N1
	
	! Writes the length of the N2 in the 2D array 
	write(1, *) N2
	
	! Writes entire array up to N1 and N2 in both dimensions
	write(1, *) Array(1:N1,1:N2)
	
	close(1)
	
	end subroutine TwoD_write
end module Arrays
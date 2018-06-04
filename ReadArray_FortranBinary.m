function A = ReadArray_FortranBinary(filename,D)
% This function reads in an D-dimensional array from unformatted FORTRAN
% Skips over record length information to properly extract data
%
% filename = 'myfile.dat'
% D = dimension of array (n=1 for 1D, n=2 for 2D ...)
%
% Example input file format for 2D array:
%
% Line | Entry         |    Data Type
%------------------------------------
%   1  | nrows         |     int32 
%   2  | ncols         |     int32 
%   3  | A(1,1)        |     double 
%   4  | A(2,1)        |     double
%   .  |   .           |       .
%   .  |   .           |       .
%   .  |   .           |       .
%  end | A(nrows,ncols)|     double

% Initializes sz row vector up to dimension length (used to store the size
%of each dimension later
sz = zeros(1,D);

% Opens the filename for binary read access, stores ID in
fileID = fopen(filename,'rb');

% Moves 4 bytes (32 binary digits) from the origin and starts at that 
% position in the data file (Fortran's data system adds 4 bytes of header
% information and 4 bytes of cloesr infromation inbetween each write
fseek(fileID, 4, 'cof');

%Loops through for each dimension, assigning the length of each dimension,
% and putting it into the sz vector
for i = 1:D
    
    % Read in length of array's ith dimension, and assigns the length to 
    % the column vector with 1 element
    sz(i) = fread(fileID, 1, 'int32');
    
    % Skips the last 4 bytes of the ith array and the first 4 bytes of
    % the ith+1 array, then starts at that new position
    fseek(fileID, 8, 'cof');
end

% Assigns A as a 1D column vector containing data
A = fread(fileID,inf,'double');

% Closes the opened file
fclose(fileID);

% Reshapes A into an array with dimensions sz(1), sz(2), up to sz(ith)
% dimension
A = reshape(A,sz);

% Transposes the matrix, since Fortran reads by columns instead of rows.
A = A';

clear all;
clc;

% Reads in image
ImageA = imread('MSU_02.jpg');
%imshow(ImageA)

% Grabs size of image
size(ImageA);
% Red = ImageA(:,:,2)

% Converts to grayscale
Inte = 0.2989 * ImageA(:,:,1) + 0.5870 * ImageA(:,:,2) + 0.1140 * ImageA(:,:,3);

% Plots current img
% surf(Inte);

% Turns to 1 inside MSU, 0 outside MSU
A = (1-Inte/255);
% surf(A);


%fid = fopen('MSU.dat', 'w');
%fwrite(fid, A);
%close(fid);

% Diffusion Method
%%change in time
dt = 4e-3;
%initial time
t0 = 0.1;
%final time
tf = .18;
%equal spacing interval for time
time = t0:dt:tf;

%Spacing
dx = 0.1;
dy = 0.1;
Nx = 60;
Ny = 118;

% Initializes matrix (allocates)
Astore = zeros(length(time),Nx,Ny);
cstore = zeros(length(time),Nx,Ny);

%Central Difference Formula for i and j directions
central_diff_i = @ (A) (circshift(A,[0 -1]) - 2.*A + circshift(A,[0 1]))./dx.^2; 
central_diff_j = @ (A) (circshift(A,[1 0]) - 2.*A + circshift(A,[-1 0]))./dy.^2; 

for n = 1:length(time)
    %new concentration gradiennt is found through old concetration value and 
    %finite difference discretization. Circshift is used for periodic boundary
    %conditions
    A = A + dt*(central_diff_i(A)+central_diff_j(A));

    %After each new concentration gradient is found, the vector is stored
    %to be graphed in real time later
    Astore(n,:,:) = A;
end

% % Plots the boundary layer smoothing
% for i = 1:length(time)
%     rshape = reshape(Astore(i,:,:),[Nx Ny]);
%     pcolor(rshape);
%     axis([0 Ny 0 Nx]);
%     pause(0.2);
% end


% Initializes concentration matrix
c = zeros(Nx,Ny);

epsilon = 1E-6;
% Makes left side all 1's, rest still zeros
c(:,1) = 1;

% So there isn't divided by 0 later on
A = A + epsilon;

% Converts all matrices to same precision
A = double(A);
Astore = double(Astore);
c = double(c);

% SBM
for n = 1:length(time)
    Areshape = reshape(Astore(n,:,:),[Nx Ny]);
    c = c + (dt./Areshape).*((Ficks(c,Areshape,0,dx,dy)-Ficks(c,Areshape,1,dx,dy))/dx+(Ficks(c,Areshape,2,dx,dy)-Ficks(c,Areshape,3,dx,dy))/dy);
    cstore(n,:,:) = c;
end

% Plots
for i = 1:length(time)
    rshape = reshape(cstore(i,:,:),[Nx Ny]);
    pcolor(rshape)
    axis([0 Ny 0 Nx]);
    pause(0.1);
end

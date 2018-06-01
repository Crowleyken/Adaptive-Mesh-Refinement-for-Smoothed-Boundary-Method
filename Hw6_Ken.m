%% CahnHilliard 1D

clc;
clear;

%Givens
W = 1;
dx = 0.1;
M = 1;
tf = 50;
e = 0.1;
N = 200;
t0 = 0;
dt = 0.001;

%Equal dts
time = t0:dt:tf;

%Allocates memory for vectors
cstore = zeros(length(time),N);
f = zeros(length(time),N);

%Initial C, random distribution of between .4 and .6, vector up to N indices
c = 0.5 + (.1)*(rand(1,N)-0.5)*2;

%Central Difference Formula
central_diff = @ (c) (circshift(c,-1) - 2.*c + circshift(c,1))./dx.^2; 

dfdc = @ (c) (W./2).*c.*(1-c).*(1-2.*c);
func = @ (c) dfdc(c) - e.^2.*central_diff(c);


for n = 1:length(time)
    c = c + M*dt*central_diff(func(c));
    cstore(n,:) = c;
    f(n,:) = (W./4).*c.^2.*(1-c).^2;
end

for i = 1:length(time) 
    plot(cstore(i,:),f(i,:),'b.');
    hold on;
    axis([-.4 1.6 0 .04]);
    hold off;
    pause(0.01);
end

%%
%Cahn-Hilliard 2D

clc
clear

%Parameters
W = 1;
M = 1;
e = 0.1;

%Grid Size
Nx = 200;
Ny = 200;

%Time restriction
t0 = 0;
tf = 5;
dt = 0.0003;
time = t0:dt:tf;

dx = 0.1;
dy = 0.1;

%x and y equal grid points
x = 0:dx:Nx;
y = 0:dy:Ny;

%allocates memory for 2D storing array
cstore = zeros(length(time),Nx,Ny);

%C is 2x2 concentration gradient where C_ave is based around 0.5
c = 0.5 + 0.1*(rand([Nx Ny])-0.5)*2;

%Central Difference Formula for i and j directions
central_diff_i = @ (c) (circshift(c,[0 -1]) - 2.*c + circshift(c,[0 1]))./dx.^2; 
central_diff_j = @ (c) (circshift(c,[1 0]) - 2.*c + circshift(c,[-1 0]))./dy.^2; 

%Functions
dfdc = @ (c) (W./2).*c.*(1-c).*(1-2.*c);
func = @ (c) dfdc(c) - e.^2.*(central_diff_i(c) + central_diff_j(c));

for n = 1:length(time)
    c = c + M*dt*(central_diff_i(func(c))+central_diff_j(func(c)));
    cstore(n,:,:) = c;
end

%allocates memory for storing array
rshape = zeros(Nx,Ny);

% rshpetest = reshape(cstore(1,:,:),[100 100])
% pcolor(rshpetest)

for i = 1:100:length(time)
    rshape = reshape(cstore(i,:,:),[Nx Ny]);
    pcolor(rshape)
    hold on;
    colormap jet;
    shading flat;
    axis equal;
    hold off;
    pause(0.001);
end
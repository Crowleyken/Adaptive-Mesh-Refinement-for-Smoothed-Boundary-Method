%Kendell Crowley
%Lab 3 Question 1, PDF of the Lab can be found in README File.

%ensuring previously stored data in MATLAB workspace is cleared
clc;
clear;

%Length
L = 20;
%Number of grids
N = 201;
%Equal spacing of x
x = linspace(0,L,N)

length(x)
%%change in x
dx = L/(N-1);
%%change in time
dt = 4e-3;
%initial time
t0 = 0.1;
%final time
tf = 1;
%initial x (center of rod)
x0 = 10;
%equal spacing interval for time
time = t0:dt:tf;

length(time)

%anonymous function for the analytical concentration gradient, 
%in terms of x,t,t0,x0, where x and t are variables (inputs a vector for x)
a_Conc = @(x,t) sqrt(t0/t).*exp(-(x-x0).^2/(4*t));

%initial concentration vector found to be used before iteration
c = a_Conc(x,t0);

%Iteration which finds the new concentration for each time step, and stores
%the vector into a 2D matrix for future use
for n = 1:length(time)
    
    %new concentration gradiennt is found through old concetration value and 
    %finite difference discretization. Circshift is used for periodic boundary
    %conditions
    c = c + dt*(circshift(c,-1) - 2*c + circshift(c,1))/dx^2;
    
    %After each new concentration gradient is found, the vector is stored
    %to be graphed in real time later
    cstore(n,:) = c;
end

%plots both the analytical and DiffEq concetration vs x over a period of
%time
for i = 1:length(time) 
    %plots the concentration gradient calculated through finite difference/
    %DiffEq
    plot(x,cstore(i,:),'b.');
    hold on;
    
    %plots the analytical concentration as it varies with time to compare
    %with DiffEq method above
    plot(x,a_Conc(x,time(i)),'r-');
    hold off;
    axis([0 L 0 1]);
    pause(0.01);
end








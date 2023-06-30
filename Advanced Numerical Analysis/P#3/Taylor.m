%% dy/dx = f    &   h = step of the grid  & i = iteration
%% y = deformation of beam   &   x = length of the beam
%% d^2y/dx^2 = g
clear;clc;
h=input ('Enter the Step Size : ');
%Boundray Conditions
x(1)=0;
e(1)=0;
y(1)=0;
%Defining the d(etta)/dx
f=@(x,e,i)-(((x(i)^4)/41750000000-(17*x(i)^3)/10437500000+(x(i)^2)/417500000+1/20875)...
    *(1+e(i)^2)^1.5);
%Definig the d^2(etta)/d x^2
g=@(x,e,i) -((9.58*10^(-11)*x(i)^3-4.887*10^(-9)*x(i)^2+4.79*10^(-9)*x(i))*(1+e(i)^2)^(3/2)+ ...
    ((2.395*10^(-11)*x(i)^4-1.629*10^(-9)*x(i)^3+2.395*10^(-9)*x(i)^2+4.79*10^(-5))^2)*(3*e(i)*(1+e(i)^2)^2));
for i=1:1:50/h
    x(i+1)=i*h;
    %Using Taylor Series Method
    e(i+1)=e(i)+h*f(x,e,i)+((h^2)/2)*g(x,e,i);
    y(i+1)=y(i)+h*e(i)+((h^2)/2)*f(x,e,i);
end
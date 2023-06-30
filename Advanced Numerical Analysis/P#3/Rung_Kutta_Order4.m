%% dy/dx = etta    &   h = step of the grid  & i = iteration
%% y = deformation of beam   &   x = length of the beam
%% ep = etta predict  &  yp = y predict
clc;clear;
h=input ('Enter the Step Size : ');
f=@(x,e,i)-((x(i)^4/41750000000-(17*x(i)^3)...
    /10437500000+x(i)^2/417500000+1/20875)*(1+e(i)^2)^1.5);
%Boundray Conditions
x(1)=0;
e(1)=0;
y(1)=0;
for i=1:1:50/h
    x(i+1)=i*h;
    k1=h*e(i);
    l1=h*f(x,e,i);
    x2(i)=x(i)+0.5*h;
    e2(i)=e(i)+0.5*l1;
    k2=h*e2(i);
    l2=h*f(x2,e2,i);
    x3(i)=x(i)+0.5*h;
    e3(i)=e(i)+0.5*l2;
    k3=h*e3(i);
    l3=h*f(x3,e3,i);
    x4(i)=x(i)+h;
    e4(i)=e(i)+l3;
    k4=h*e4(i);
    l4=h*f(x4,e4,i);
    e(i+1)=e(i)+(l1+2*l2+2*l3+l4)/6;
    y(i+1)=y(i)+(k1+2*k2+2*k3+k4)/6;
end
    
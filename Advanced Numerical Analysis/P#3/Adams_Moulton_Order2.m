%% dy/dx = etta    &   h = step of the grid  & i = iteration
%% y = deformation of beam   &   x = length of the beam
%% ep = etta predict  &  yp = y predict
clc;clear;
h=input ('Enter the Step Size : ');
f=@(x,e,i)(x(i)^4/41750000000-(17*x(i)^3)/10437500000+x(i)^2/417500000+1/20875)*(1+e(i)^2)^1.5;
fp=@(x,ep,i)(x(i)^4/41750000000-(17*x(i)^3)/10437500000+x(i)^2/417500000+1/20875)*(1+ep(i)^2)^1.5;
%Boundray Conditions
x(1)=0;
e(1)=0;
y(1)=0;
for i=1:1:50/h
    x(i+1)=(i)*h;
    %Using Modified Euler to Predict the first Points
    yp(i+1)=y(i)+h*e(i);
    ep(i+1)=e(i)+h*f(x,e,i);
    y(i+1)=y(i)+(h/2)*(e(i)+ep(i+1));
    e(i+1)=e(i)+(h/2)*(f(x,e,i)+fp(x,ep,i+1));
    %Using Adams Moulton Method Second Order 
    if i>=3
        yp(i+1)=y(i)+(h/12)*(23*e(i)-16*e(i-1)+5*e(i-2));
        ep(i+1)=e(i)+(h/12)*(23*f(x,e,i)-16*f(x,e,i)+5*f(x,e,i-2));
        y(i+1)=y(i)+(h/12)*(5*ep(i+1)+8*e(i)-e(i-1));
        e(i+1)=e(i)+(h/12)*(5*fp(x,ep,i+1)+8*f(x,e,i)-f(x,e,i-1));
    end
end
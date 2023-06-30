%% u:Wavelength
%% T:Temprature of blackbody                                                     
%% x:The emissive power of Blackbody per wavelength 
tic
clc 
clear
syms x h;
format long

%Temprature of Blackbody
T=1000;

%Defining amount of Emissive Power
u=100;

%Definig a 3 variable Function
f=@(T,u,x)u-((3.74e+08)/((x^5)*(exp(14400/(x*T))-1)));

% Epsilon
e=0.000001;

% Initial guess for emissive power
x_int=[1,0.5];

%Number of Iterations
k=1;

x0=x_int(1);
m=@(x)(f(T,u,x));
%find the derivative
g=inline(limit(((m(x+h)-m(x))/h),h,0));
%Using first degree Newton Method
xn=x0-(m(x0)/g(x0));
%Stop the Iterations
    while abs(xn-x0)>e
      %Find & Replace new Solutions
       x0=xn;
       xn=x0-(m(x0)/g(x0));
       k=k+1;
    end
%Results
disp (['Number of iterations are : ',num2str(k)]);
disp (['Suppose that Emissive power is ',num2str(u),' and Temprature is ',num2str(T),' K']);
fprintf ('%s','The amount of Wavelength is : ');fprintf('%d\n',xn);
toc
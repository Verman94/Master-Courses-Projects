%% u:Wavelength
%% T:Temprature of blackbody                                                     
%% x:The emissive power of Blackbody per wavelength 
tic
clc 
clear
format long

%Temprature of Blackbody
T=1000;

%Defining amount of emissive power
u=100;

%Definig a 3 variable Function
f=@(T,u,x)u-((3.74e+08)/((x^5)*(exp(14400/(x*T))-1)));

% Epsilon
e=0.000001; 

% Initial guess for emissive power
x_int=[1,0.5];

%number of Iterations
 k=1;  
 
x0=x_int(1);
x1=x_int(2);
m=@(x)(f(T,u,x));

%using Bisection Method
if m(x0)*m(x1)<0
    xn=(x0+x1)/2;
else
    disp 'Your initial guess is wrong!'
    break;
end
while abs(xn-x0)>e
    k=k+1;
    if m(xn)*m(x0)<0
        x1=xn;
        xn=(x1+x0)/2;
    elseif m(xn)*m(x1)<0
        x0=xn;
        xn=(x0+x1)/2;
    end
end
%Results
disp (['Number of iterations are :  ',num2str(k)]);
disp (['Suppose that Emissive power is ',num2str(u),' and Temprature is ',num2str(T),' K']);
fprintf ('%s','The amount of Wavelength is : ');fprintf('%d\n',xn);
toc

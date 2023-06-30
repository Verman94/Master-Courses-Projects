%% u:Wavelength
%% T:Temprature of blackbody                                                     
%% x:The emissive power of Blackbody per wavelength 
tic
clc 
clear
format long
syms x h;

%Temprature of Blackbody
T=1000;

%Defining amount of Emissive Power
u=100;


%Definig a 3 variable Function
%g=@(T,u,x)(x^5*exp(14400/x/T)-u/3.74e+08)^(1/5);
g=@(T,u,x)14400/(log((3.74e+08/(x^5*u))*T));
%g=@(T,u,x)(x^5*u*exp(14400/(x*T))-3.74e+08)/(x^4*u);
%g=@(T,u,x)3.74e+08/(x^4*u*(exp(14400/(u*T))-1));

% Epsilon
e=0.000001; 

%Number of iterations
k=1;

% Initial guess for emissive power
x_int=[4,2];
x0=x_int(1);
m=@(x)(g(T,u,x));
f=inline(limit(((m(x+h)-m(x))/h),h,0));

%Using Fixed Point Method
xn=m(x0);
for i=1:k
    k=k+1;
while abs(xn-x0)>e
      x0=xn;
      xn=m(x0);
      landa(1,i)=xn;
      if k>500
          break;
      end
end
end
%Results
disp (['Number of iterations are :  ',num2str(k)]);
disp (['Suppose that Emissive power is ',num2str(u),' and Temprature is ',num2str(T),' K']);
fprintf ('%s','The amount of Wavelength is : ');fprintf('%f\n',xn);
toc
%% u:Wavelength
%% T:Temprature of blackbody                                                     
%% x:The emissive power of Blackbody per wavelength
tic
clc 
clear
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
x_int=[1,2,3]; 
x2=x_int(3);
x1=x_int(2);
x0=x_int(1);
h1=x1-x0;
h2=x0-x2;
Y=h2/h1;
m=@(x)(f(T,u,x));

%Using the Muller Method
a=(m(x1)-Y*m(x0)+m(x2))/(Y*(h1)^2);
b=(m(x1)-m(x0)-(a*(h1)^2))/h1;
c=m(x0);

%Deciding the sign of b+-sqrt(b^2-4ac)
if b>=0
  xn=x0-(c/(b^2+sqrt(b^2-4*a*c)));
else
  xn=x0-(c/(b^2-sqrt(b^2-4*a*c)));
end
%Number of Iterations
k=1;
%Stop the Iteration
while abs(xn-x0)>=e
k=k+1;
%Determine the next three points
    if xn>x0
      x2=x0;
      x0=xn;
    else
       x1=x0;
       x0=xn;
    end
   %Find & Replace new Solutions
    h1=(x1-x0);
    h2=x0-x2;
    Y=h2/h1;
    a=(m(x1)-Y*m(x0)+m(x2))/(Y*(h1)^2);
    b=(m(x1)-m(x0)-(a*(h1)^2))/h1;
    c=m(x0); 
    if b>=0
      xn=x0-(c/(b^2+sqrt(b^2-4*a*c)));
    else
      xn=x0-(c/(b^2+sqrt(b^2-4*a*c)));
    end 
end
disp (['Number of iterations are :  ',num2str(k)]);
disp (['Suppose that Emissive power is ',num2str(u),' and Temprature is ',num2str(T),' K']);
fprintf ('%s','The amount of Wavelength is : ');fprintf('%d\n',xn);
toc
%Solving 2D Heat Transfer Equation for SHAPE 1 Using Implicit Method
clear all
clc
%Used Parameters
a = 4e-2;        %  Lenght of SHAPE1  [m]
b = 4e-2;        %  Width of SHAPE1  [m]
T_inf1 = 25;     %  Temprature of upper outside  [C]
T_inf2 = 20;     %  Temprature of Underside outside  [C]
q = 25;          %  Rate of Generate Heat [W]
h_conv1 = 45;    %  Heat Convection Coefficient of upper Air  [W/m^2.K]
h_conv2 = 35;    %  Heat Convection Coefficient of underside Air  [W/m^2.K]
k = 0.6;  %  Heat Conduction Coefficient of Water  [W/m.C]
alfa = 0.143e-6; %  Thermal Diffusivity of Water [m^2/s]
t = 3600;           %  Time [s]
%Numerical Parameters
dx = 0.002;
dy = 0.002;
nx = round(a / dx + 1);  %  size of the X-Grid
ny = round(b / dy + 1);  %  size of the Y-Grid
% define the time step
dt = 0.4;
l = round(t / dt + 1);  %  number of the iteration
Sx = (alfa*dt)/dx^2;  %  Coefficient in Heat equation
Sy = (alfa*dt)/dy^2;  %  Coefficient in Heat equation
C1=2*h_conv1*dy/k;C2=2*h_conv2*dy/k;
C3=(q*alfa*dt)/k;
%Set the Dirichlet Boundray Condition
IT(1,1:ny) = 100;
%Set the Temprature at initial time
IT(2:nx,1:ny) = 50;
%temp=cell(1,l);
T=IT';
%temp{1}=T;         %Using this cell array for diffrent times
A=zeros((nx-1)*ny,(nx-1)*ny);   % Because one side has dirichlet B.C
tic
for p=2:l
for i=1:(nx-1)*ny
    A(i,i)= -2*Sx-1-2*Sy;
end
for i=1:ny*(nx-2)
    A(i,i+ny)=Sx;
end
for i=ny*(nx-2)+1:ny*(nx-1)
    A(i,i-ny)=2*Sx;
end
for i=1:ny*(nx-3)
    A(i+ny,i)=Sx;
end
for j=1:nx-1
for i=(j-1)*ny+2:j*ny-1
    A(i,i-1)=Sy;
    A(i,i+1)=Sy;
end
end
for i=1:nx-1
    A((i-1)*ny+1,(i-1)*ny+2)=Sy*((1/(C1-1))+1);      %Top Boundray Layer
    A(i*ny,i*ny-1)=Sy*(1/(C2+1)+1);                  %Bottom Boundray Layer
end
%Matrix of Constants
for i=1:nx-1
    for j=1:ny
        if i==1
            if j==1
                B(i,j)=-IT(i+1,j)-C3-Sx*IT(i,j)+Sy*T_inf1*C1/(C1-1);  %Top Left Edge
            elseif j==ny
                B(i,j)=-IT(i+1,j)-C3-Sx*IT(i,j)-Sy*C2*T_inf2/(C2+1);  %Bottom left edge
            else
                B(i,j)=-IT(i+1,j)-C3-Sx*IT(i,j); %first row interior nodes
            end
        else
            if j==1
                B(i,j)=-IT(i+1,j)-C3+Sy*T_inf1*C1/(C1-1);  %Top Edge
            elseif j==ny
                B(i,j)=-IT(i+1,j)-C3-Sy*C2*T_inf2/(C2+1);  %Bottom Edge
            else
                B(i,j)=-IT(i+1,j)-C3;  %Interior nodes
            end
        end
    end
end
b=reshape(B',[],1);
Q=Seidel(A,b,0.1);   %Using gauss seidel with relaxation factor
IT(1,1:ny)=100;
for i=2:nx
    for j=1:ny
        IT(i,j)=Q(j+(i-2)*ny,1);
    end
end
mnm=IT';
%temp{p}=IT';
end
for i=1:nx
    X((i-1)*ny+1:i*ny,1)=(i-1)*dx;
end
for i=1:ny
    for j=1:nx
        Y(((i-1)*ny)+j,1)=(j-1)*dy;
    end
end
T=reshape(mnm,[],1);
toc
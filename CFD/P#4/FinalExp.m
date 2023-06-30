%Solving 2D Heat Transfer Equation for SHAPE 1 Using Explicit Method
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
k = 0.6;         %  Heat Conduction Coefficient of Water  [W/m.C]
alfa = 0.143e-6; %  Thermal Diffusivity of Water [m^2/s]
t = 3600;           %  Time [s]
%Numerical Parameters
dx = input ('Enter the Step in direction of X in meter = ');
dy = input ('Enter the Step in direction of Y in meter = ');
nx = round(a / dx + 1);  %  size of the X-Grid
ny = round(b / dy + 1);  %  size of the Y-Grid
% define the time step
%dt = (dx^2+dy^2) / alfa / 8;
dt=0.4;
l = round(t / dt + 1);  %  number of the iteration
Sx = (alfa*dt)/dx^2;  %  Coefficient in Heat equation
Sy = (alfa*dt)/dy^2;  %  Coefficient in Heat equation
C1=2*h_conv1*dy/k;C2=2*h_conv2*dy/k;
%Set the Dirichlet Boundray Condition
IT(1,1:ny) = 100;
%Set the Temprature at initial time
IT(2:nx,1:ny) = 50;
temp=cell(1,l);
temp{1}=IT';         %Using this cell array for diffrent times
% Using Explicit Method Formulation to find temprature in each time step
for p=2:l
    %Middle Points of the Grid
    for i=2:nx-1
        for j=2:ny-1
            Q(i,j)=IT(i,j)+Sx*(IT(i+1,j)-2*IT(i,j)+IT(i-1,j))+Sy*(IT(i,j+1)-...
                2*IT(i,j)+IT(i,j-1))+(alfa*dt*q)/k;
        end
    end
    %Top Edge of the Grid
    for i=2:nx
        j=1;
        if i~=nx
            Q(i,j)=IT(i,j)+(q*alfa*dt/k)+Sx*(IT(i+1,j)-2*IT(i,j)+IT(i-1,j))+...
                Sy*((1/(C1+1))*(IT(i,j+1)+C1*T_inf1)-2*IT(i,j)+IT(i,j+1));
        else 
            %Top Right Corner of the Grid
            Q(i,j)=IT(i,j)+(q*alfa*dt/k)+Sx*(2*IT(i-1,j)-2*IT(i,j))+...
                Sy*((1/(C1+1))*(IT(i,j+1)+C1*T_inf1)-2*IT(i,j)+IT(i,j+1));
        end
    end
    %Bottom Edge of the Grid
    for i=2:nx
        j=ny;
        if i~=nx
            Q(i,j)=IT(i,j)+(q*alfa*dt/k)+Sx*(IT(i-1,j)-2*IT(i,j)+IT(i+1,j))+...
                Sy*((1/(C2+1))*(C2*T_inf2+IT(i,j-1))-2*IT(i,j)+IT(i,j-1));
        else
            %Bottom Right Corner of the Grid
            Q(i,j)=IT(i,j)+(q*alfa*dt/k)+Sx*(2*IT(i-1,j)-2*IT(i,j))+...
                Sy*((1/(C2+1))*(C2*T_inf2+IT(i,j-1))-2*IT(i,j)+IT(i,j-1));
        end
    end
    %Right side of the Grid
    for j=2:ny-1
        i=nx;
            Q(i,j)=IT(i,j)+(q*alfa*dt/k)+Sx*(2*IT(i-1,j)-2*IT(i,j))+...
                Sy*(IT(i,j+1)-2*IT(i,j)+IT(i,j-1));
    end
    %Left Side of the Grid
    Q(1,1:ny)=100;
    IT=Q;
    temp{p}= IT';
end
for i=1:nx
    X((i-1)*ny+1:i*ny,1)=(i-1)*dx;
end
for i=1:ny
    for j=1:nx
        Y(((i-1)*ny)+j,1)=(j-1)*dy;
    end
end
T=reshape(temp{l},[],1);
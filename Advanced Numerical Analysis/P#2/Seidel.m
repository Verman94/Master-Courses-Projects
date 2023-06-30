%% A = Matrix of Coefficients
%% B = Matrix of Constants
%% eps = amount of Tolerance to Stop the Iteration
%% k = Number of Iterations
%% alfa = Relaxation Factor

clear all
clc
%Asking the User to Give us the Matrix
A=[-12	1	0	0	0	0	0	0	1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0
1	-24	1	0	0	0	0	0	0	2	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0
0	1	-24	1	0	0	0	0	0	0	2	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0
0	0	1	-24	1	0	0	0	0	0	0	2	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0
0	0	0	1	-24	1	0	0	0	0	0	0	2	0	0	0	0	0	0	0	0	0	0	0	0	0	0
0	0	0	0	1	-24	1	0	0	0	0	0	0	2	0	0	0	0	0	0	0	0	0	0	0	0	0
0	0	0	0	0	1	-24	1	0	0	0	0	0	0	2	0	0	0	0	0	0	0	0	0	0	0	0
0	0	0	0	0	0	1	-12	0	0	0	0	0	0	0	1	0	0	0	0	0	0	0	0	0	0	0
1	0	0	0	0	0	0	0	-4	2	0	0	0	0	0	0	1	0	0	0	0	0	0	0	0	0	0
0	1	0	0	0	0	0	0	1	-4	1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0
0	0	1	0	0	0	0	0	0	1	-4	1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0
0	0	0	1	0	0	0	0	0	0	1	-4	1	0	0	0	0	0	0	0	0	0	0	0	0	0	0
0	0	0	0	1	0	0	0	0	0	0	1	-4	1	0	0	0	0	0	0	0	0	0	0	0	0	0
0	0	0	0	0	1	0	0	0	0	0	0	1	-4	1	0	0	0	0	0	0	0	0	0	0	0	0
0	0	0	0	0	0	1	0	0	0	0	0	0	1	-4	1	0	0	0	0	0	0	0	0	0	0	0
0	0	0	0	0	0	0	1	0	0	0	0	0	0	2	-4	0	0	0	0	0	0	0	0	0	0	0
0	0	0	0	0	0	0	0	1	0	0	0	0	0	0	0	-4	2	0	0	0	0	0	0	1	0	0
0	0	0	0	0	0	0	0	0	1	0	0	0	0	0	0	1	-4	1	0	0	0	0	0	0	1	0
0	0	0	0	0	0	0	0	0	0	2	0	0	0	0	0	0	2	-10	1	0	0	0	0	0	0	1
0	0	0	0	0	0	0	0	0	0	0	2	0	0	0	0	0	0	1	-8	0	0	0	0	0	0	0
0	0	0	0	0	0	0	0	0	0	0	0	2	0	0	0	0	0	0	1	-8	0	0	0	0	0	0
0	0	0	0	0	0	0	0	0	0	0	0	0	2	0	0	0	0	0	0	1	-8	0	0	0	0	0
0	0	0	0	0	0	0	0	0	0	0	0	0	0	2	0	0	0	0	0	0	1	-8	0	0	0	0
0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	1	0	0	0	0	0	0	1	-4	0	0	0
0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	1	0	0	0	0	0	0	0	-2	1	0
0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	2	0	0	0	0	0	0	1	-4	1
0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	1	0	0	0	0	0	0	1	-4
];
b=[-7000
-14000
-14000
-14000
-14000
-14000
-14000
-68000
0
0
0
0
0
0
0
0
0
0
-1200
-1200
-1200
-1200
-1200
-600
0
0
-600
];
% A=input('Enter your Coefficients Matrix :');
% b=input('Enter your Constant Matrix :');
DA = length(A);
% Initial Value
x=zeros(DA,1);
% Tolerance Value
epsilon=0.00001;   
k=1;
tole=1;
%Relaxation Factor
alpha=input('Enter your Relaxation Factor :');
%Using Gauss Seidel Algorithem
while tole>epsilon  
for i=1:DA
    sigma1=0;
    if i~=1
    for m=1:i-1
       sigma1=sigma1+A(i,m)*x(m,k+1);
    end
    end
    sigma2=0;
    for j=i:DA
        sigma2=sigma2+A(i,j)*x(j,k);
    end   
    x(i,k+1)=x(i,k)+(alpha)*(b(i)-sigma1-sigma2)/A(i,i);
    R(i,1)=abs((b(i)-sigma2-sigma1)/A(i,i));
end
%Number of iterations
 k=k+1;
 %Calculate the sum of Residuals
 sum = 0;
 for j=1:DA
     sum=sum+abs(R(j));
 end
tole=sum/DA;
end
X=x(1:DA,k);
disp(['X= ',mat2str(X)]);
disp(['Number of iterations : ',num2str(k)]);
%% A=Coefficient Matrix
%% B=Constant Matrix
%% DA=Columns of Matrix A
%% DB=Rows of Matrix B
%% AX=B  &  X=(INVERSE_A)*B  &  A*(INVERSE_A)=I
%% U=Upper Triangular Matrix 
%% L=Lower Triangular Matrix
%% A=LU  &  U*(INVERSE_A)=C  &  LC=I

clear all
clc
%Asking the User to Give us the Matrix
A=input ('insert the Matrix of Coefficients : ');
B=input ('insert the Matrix of Constants : ');
%to findout how many equations we have
DA=size(A,2);DB=size(B,1);
%Investigate if the Matrixes are in same dimension
if DA~=DB ,disp('Dimension of your Entered Matrixes are not similar!');end
for m=1:DA-1
    %investigate if Matrix needs pivoting
    [max_in_column,i]=max(abs(A(m:DA,m)));
    ipr=i+m-1;
        if ipr~=m
            %pivoting
            A([m,ipr],:)=A([ipr,m],:);
            B([m,ipr],:)=B([ipr,m],:);
        end
%Find L & U
U=zeros(DA);L=zeros(DA);
L(:,1)=A(:,1);
U(1,:)=A(1,:)/L(1,1);
U(1,1)=1;
for k=2:DA  
    for j=2:DA
        for i=j:DA
            L(i,j)=A(i,j)-dot(L(i,1:j-1),U(1:j-1,j));
        end
        U(k,j)=(A(k,j)-dot(L(k,1:k-1),U(1:k-1,j)))/L(k,k);
    end
end
%Creating a Unit Matrix
I=eye(DA);
%Using L-U Decomposition to find inverse A
for z=1:DA
    a=I(:,z);
    c(1)=a(1)/L(1,1);
    for i=2:DA
        S=0;
        for k=1:i-1
            S=S+L(i,k)*c(k);
        end
        S=a(i)-S;
        c(i)=S/L(i,i);
    end
    %Find the X 
    X=zeros(DA,1);
    for i=DA:-1:1
        S=0;
        for j=i+1:DA
            S=S+U(i,j)*X(j);
        end
        X(i)=c(i)-S;
    end
    h(z,:)=X;
    Inverse_A=h';
end
end
%RESULTS
X=Inverse_A*B;
disp('X = ');disp(mat2str(X));
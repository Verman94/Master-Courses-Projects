%% A=Coefficient Matrix
%% B=Constant Matrix
%% DA=Columns of Matrix A
%% DB=Rows of Matrix B
%% AX=B

clear all
clc
%Asking the User to Give us the Matrix
A=input ('insert the Matrix of Coefficients : ');
B=input ('insert the Matrix of Constants : ');
%to findout how many equations we have
DA=size(A,2);DB=size(B,1);
%Investigate if the Matrixes are in same dimension
if DA~=DB ,disp('Dimension of your Entered Matrixes are not similar!');end
for k=1:DA-1
    %investigate if Matrix is Singular and pivoting
    [max_in_column,i]=max(abs(A(k:DA,k)));
    ipr=i+k-1;
        if ipr~=k
            %changing rows
            A([k,ipr],:)=A([ipr,k],:);
            B([k,ipr],:)=B([ipr,k],:);
        end
    %Use Gaussian Elimination Algorithem
    for i=k+1:DA
        b(k)=A(i,k);
        c(k)=A(k,k);
        e(k)=B(i);
        f(k)=B(k);
        e(k+1)=e(k)-(b(k)./c(k)).*f(k);
        B(i)=e(k+1);
        for j=k+1:DA
            a(k)=A(i,j);
            d(k)=A(k,j);
            a(k+1)=a(k)-(b(k)./c(k)).*d(k);
            A(i,j)=a(k+1);
        end
    end
end
A=triu(A);
%Find the amount of X
X(DA)=B(DA)/A(DA,DA);
for i=DA-1:-1:1
    c=0;
    for j=i+1:DA
        c=c+(A(i,j)*X(j));
    end
    X(i)=(B(i)-c)/A(i,i);
end
%RESULT
disp('X = ');disp(mat2str(X));
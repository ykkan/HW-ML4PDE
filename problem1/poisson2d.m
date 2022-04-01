function U = poisson2d(f, n)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
K1D = spdiags(ones(n,1)*[-1 2 -1],-1:1,n,n);
I1D = speye(size(K1D));
A = kron(K1D,I1D)+kron(I1D,K1D);
x=linspace(0,1,n);
y=x;
h=x(2)-x(1);
rhs=zeros(n^2,1);
for j=1:n
    for i=1:n
        k = i+(j-1)*n;
        rhs(k)=f(x(i),y(j))*h^2;
    end
end
U=A\rhs;
end
function [L,U]=geMLU(N,u,v)
%
% This function computed A=L*U by GTH-like algorithm of an M-matrix A given by
% a triplet representation (N, u, v), where
%
%  N  such that A(i,j)=-N(i,j) for i .ne. j.; N(i,j) not used.
%  A*u=v. u positive vector, and v nonnegative vector
%
% GTH-like algorithm is given in
%
%   Attahiru Sule Alfa, Jungong Xue and Qiang Ye.
%   Accurate computation of the smallest eigenvalue of a diagonally dominant $M$-matrix.
%   Math. Comp. 71 (2002) 217-236.
%
%             RCL 10/19/2019
%---------------------------------------------------
%
% Input
%
%       N     (n-by-n) matrix
%             A(i,j)=-N(i,j) for i .ne. j
%
%       u     (n-by-1) positive vector
%
%       v     (n-by-1) nonegative vector
%             A*u=v
%
% Output
%
%      L      (n-by-n) lower triangular matrix, all diagonal entries are  1 
%      U      (n-by-n) upper triangular matrix
%----------------------------------------------------

n=size(N,1);

for i=1:n-2,
    N(i,i)=( v(i)+N(i,i+1:n)*u(i+1:n) )/ u(i);
    tmp=1.0/N(i,i);
    N(i+1:n,i)=tmp*N(i+1:n,i);
    N(i+1:n,i+1:n)=N(i+1:n,i+1:n)+N(i+1:n,i)*N(i,i+1:n); % diagonal entries not useful
    v(i+1:n)=v(i+1:n)+v(i)*N(i+1:n,i);
end
%
i=n-1;
N(i,i)=( v(i)+N(i,n)*u(n) )/ u(i);
tmp=1.0/N(i,i);
N(n,i)=tmp*N(n,i);
v(n)=v(n)+v(i)*N(n,i);
N(n,n)=v(n)/u(n);

L=-tril(N)+diag(diag(N))+diag(ones(n,1));
U=-triu(N)+2*diag(diag(N));


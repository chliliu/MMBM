function X=geMx(N,u,v,B,type)
%
% This function solves A*x=b or x*A=b, by GTH-like algorithm of an M-matrix A given by
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
%             RCL 8/2/2010
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
%       B     (n-by-k) matrix if type=1; 
%             (k-by-n) matrix if type=2;
%             Although the algorithm guarantees high relative accuracy for b>=0, this code
%             does not require b>=0.
%
%       type  (integer) optional
%             default type = 1: Solve A*X=B
%                          = 2: solve X*A=B
%
% Output
%
%      X      (n-by-k) matrix if type=1; 
%             (k-by-n) matrix if type=2;
%----------------------------------------------------

n=size(N,1);
%L=diag(ones(n,1)); U=zeros(n,n);

testLU=0;  % 1: test A-L*U; 0: do not test A-L*U
if testLU==1,
   A=diag((N*u+v)./u)-N;
end

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

if testLU==1,
   % A=L*U, where
   L=-tril(N)+diag(diag(N))+diag(ones(n,1));
   U=-triu(N)+2*diag(diag(N));
   disp(norm(A-L*U,1)/norm(A,1));
end

if nargin == 4,
   type = 1;   % This is default, too.
end
if type == 1, % solve A*X=B
   % L*Y=B
   Y=zeros(size(B)); 
   Y(1,:)=B(1,:);
   for i=2:n,
       Y(i,:)=B(i,:)+N(i,1:i-1)*Y(1:i-1,:);
   end
   %
   % U*X=Y
   X=zeros(size(B));
   X(n,:)=Y(n,:)/N(n,n);
   for i=n-1:-1:1,
       X(i,:)=( Y(i,:)+N(i,i+1:n)*X(i+1:n,:) )/N(i,i);
   end
elseif type == 2, % type = 2: solve X*A=B
   % Y*U=B
   Y=zeros(size(B)); 
   Y(:,1)=B(:,1)/N(1,1);
   for j=2:n,
       Y(:,j)=( B(:,j)+Y(:,1:j-1)*N(1:j-1,j) )/N(j,j);
   end
   % X*L=Y   
   X=zeros(size(B));
   X(:,n)=Y(:,n);
   for j=n-1:-1:1,
       X(:,j)=Y(:,j)+X(:,j+1:n)*N(j+1:n,j);
   end
end

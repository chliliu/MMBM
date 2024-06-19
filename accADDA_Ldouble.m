function [X1, Y1, info]=accADDA_Ldouble(m, n, W, opts)
%
% Alternating Directional Doubling Algorithm for MARE
%
%                X D X-A X-X B + C = 0,  
% where 
%
%             m    n
%           [ B   -D ]  m
%   W =     [-C    A ]  n
%
% is a nonsingular or irredicuble singular M-matrix. As a by-product, it also solves the dual equation
%
%                 D - Y A - B Y + Y C Y = 0.
%
% This is the highly accurate version, and it is assumed a left triplet representation, in contrast to accADDA(...),
%
%    W={N,u,v} such that W(i,j)=-N(i,j)<=0 for i .ne. j, and 
%    u*W=v to recover diag(W)=(u*N+v)./u, where u is a positive row-vector
%    while v is a row-nonnegative vector.
%
% is given.
%
% An algorithm, essentially the same as 
%
% Jungong Xue and Ren-Cang Li.
% Highly accurate doubling algorithms for ${M}$-matrix algebraic {Riccati} equations.
% {\em Numer. Math.}, 135(3):733--767, 2017.
%
% can be obtained verbatim. But here we take a short cut. Namely, we solve instead
%
%   X'  D'  X' - B' X' - X' A' + C' = 0                (***)
%
% and call accADDA(...). The W-matrix of (***) is
%
%           n   m
%         [ A' -D'] n
%  W1 =   [ -C' B'] m   
%
% whose (right) triplet representation is easily obtained from the left triplet representation of W
%
%   RCL 10/19/2019
%
%---------------------------------------------------
%
% Usage:
%
%    [X1, Y1, info]=accADDA(m,n,W,opts)
%
%---------------------------------------------------
%
% Input
%
%       m     positive integer to partition W
%
%       n     positive integer to partition W
%
%       W     (structure) for the triplet representation
%             W.N: the opposite of the off-diagonal part of W. It is nonnegative                  
%                  W(i,j)=-W.N(i,j)<=0 for i .ne. j
%                  W.N(i,i)=0 for all 1<=i<=m+n
%             W.u: (m+n)-row vector, positive
%             W.v: (m+n)-row vector, nonegative
%
%       opts  (structure) 
%             opts.tol      relative entrywise tolerence
%             opts.maxitn:  number of maximal iteration steps, 
%                           default 40;
%             opts.res      =1  report residuals at each iteration
%                           =0  don't report residuals
%             opts.nrms:    [norm(A,1) norm(B,1) norm(C,'fro') norm(D,1)]
%                           needed only when opts.res = 1;
%             opts.meth     which method is used for inverting I-Xk*Yk and I-Yk*Xk
%                           = 1: GTH but with triplet representation for them simply 
%                                calculated by v_1=u_1*(I-Yk*Xk) and v_2=u_2*(I-Xk*Yk)
%                           = 2: GTH  triplet representations recursively calculated 
%                                as in Xue's & Li's paper.
%             opts.exact    = 1  exact Phi is given and ||X-Phi||/||X||_F is reported
%                           = 0  no Phi is provided
%             opts.Phi      the exact solution, provided only when opts.exact = 1.
%
% Output
%
%       X1    (n-by-m) matrix
%             minimal nonegative solution.
%
%       Y1    (m-by-n) matrix
%             minimal nonegative solution.
%
%       info  (optional) struct
%             info.itn  number of iterations taken
%             info.upd  ith-by-2 or 4
%                       (i,1:2)=[abs. ith X-update, rel. ith X-update]
%                       (i,3:4) [abs. ith Y-update, rel. ith Y-update] (if info.upd is ith-by-4)
%             info.stop quantities used to terminate the iteration (Kahan's criteria)
%                                        [ (X1-X2).^2 ]_{i,j}
%                       max ------------------------------------------------                   
%                            [ (X0-X1)_{i,j} - (X1-X2)_{i,j} ]* (X2)_{i,j}
%                        
%             info.succ 1  successful exit after satisfying the stopping tolerance
%                       0  exit after the maxitn is reached.
%             info.res  (itn+1)-by-2 
%                       only referenced when opts.res=1 for recording residuals at each iteration calculated as
%
%                       info.res(i,1) -- usual residual
%                                       norm(X *D* X-A* X-X* B + C,'fro')
%                         ---------------------------------------------------------------------------
%                             norm(X,'fro')*[ norm(D,1)*norm(X,1)+norm(A,1)+norm(B,1) ]+norm(C,'fro')
%
%                       info.res(i,2) -- Entrywise residual (see Xue & Li's paper)
%
%                       ideally, norm(...,1) should be norm(...), but it is used for computational
%                       convenience.
%             info.errX only referenced when opts.exact=1 for recording relative errors between eact approximation X and Phi
%                             norm(X-Phi,'fro')
%                         ---------------------
%                              norm(Phi,'fro')
%----------------------------------------------------

mpn=m+n; 

m1=n; n1=m;
W1.N=[W.N(m+1:mpn,m+1:mpn)' W.N(1:m,m+1:mpn)';
      W.N(m+1:mpn,1:m)' W.N(1:m,1:m)'];
W1.u=[W.u(m+1:mpn)';  W.u(1:m)'];
W1.v=[W.v(m+1:mpn)'; W.v(1:m)'];

if opts.exact==1,
   opts.Phi=opts.Phi';
end

[X1, Y1, info]=accADDAdouble(m1, n1, W1, opts);

% transpose back
if opts.exact==1,
   opts.Phi=opts.Phi';
end
X1=X1'; Y1=Y1';



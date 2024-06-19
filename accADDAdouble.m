function [X1, Y1, info]=accADDAdouble(m, n, W, opts)
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
% This is the highly accurate version, and it is assumed a triplet representation
%
%    W={N,u,v} such that W(i,j)=-N(i,j)<=0 for i .ne. j, and 
%    W*u=v to recover diag(W)=(N*u+v)./u, where u is a positive
%    vector while v is a nonnegative vector.
%
% is given.
%
% This is based on
%
% Jungong Xue and Ren-Cang Li.
% Highly accurate doubling algorithms for ${M}$-matrix algebraic {Riccati} equations.
% {\em Numer. Math.}, 135(3):733--767, 2017.
%
%   RCL 10/4/2015
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
%             W.u: (m+n)-vector, positive
%             W.v: (m+n)-vector, nonegative
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
%                                calculated by v_1=(I-Yk*Xk)*u_1 and v_2=(I-Xk*Yk)*u_2
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
%             info.eerr_X only referenced when opts.exact=1 for recording relative errors between eact approximation X and Phi
%                         max_{i,j} [(X-Phi)./Phi]_{i,j}
%
%----------------------------------------------------

mpn=m+n; 

C=W.N(m+1:mpn,1:m);       D=W.N(1:m,m+1:mpn);
N_A=W.N(m+1:mpn,m+1:mpn); N_B=W.N(1:m,1:m);
diagW=(W.N*W.u+W.v)./W.u;

u2=W.u(m+1:mpn);    u1=W.u(1:m); 
%v2=W.v(m+1:mpn);    v1=W.v(1:m);

% alpha and beta cannot be 0
alpha=1.01/max(diagW(m+1:mpn)); beta=1.01/max(diagW(1:m));

% compute initial E0, F0, X0, Y0, and w0
u0=[u1/alpha; u2/beta]; v0=W.v+u0;
Nwk=[alpha*N_B, beta*D;
     alpha*C,   beta*N_A];
Bwk=[beta*N_B, alpha*D; beta*C, alpha*N_A];
for i=1:m
    Bwk(i,i)=1.0-diagW(i)*beta;
end
for i=1:n    
    Bwk(m+i,m+i)=1.0-diagW(m+i)*alpha;
end
if opts.meth==2
   Bwk=[Bwk W.v];
end
T=geMx(Nwk,u0,v0,Bwk,1);
E0=T(1:m,1:m); Y0=T(1:m,m+1:mpn);
X0=T(m+1:mpn,1:m); F0=T(m+1:mpn,m+1:mpn);
if opts.meth==2
   w0=(alpha+beta)*T(:,mpn+1);
   w1=zeros(mpn,1);
end  

Xdiff0=X0; Ydiff0=Y0;
abserr_upd0=max(max(max(X0)),max(max(Y0)));

itn=0; 

if opts.res==1,
   nrmA=opts.nrms(1);
   nrmB=opts.nrms(2);
   nrmC=opts.nrms(3);
   nrmD=opts.nrms(4);
   
   RR=zeros(n,m); 
   for i=1:n
       RR(i,:)=diagW(m+i)*X0(i,:);
   end
   for j=1:m
       RR(:,j)=RR(:,j)+X0(:,j)*diagW(j);
   end
   RL=X0*D*X0+N_A*X0+X0*N_B+C;
   
   ResX=RL-RR;
   eres=max(max(abs(ResX)./RR)); 

   nrmX=norm(X0,'fro'); nrmXp=norm(X0,1);
   res=norm(ResX,'fro')/(nrmX*(nrmXp*nrmD+nrmA+nrmB)+nrmC);
   
   ress=[res eres];
end

if opts.exact == 1,
   nrmPhi=norm(opts.Phi,'fro'); 
   err_X=norm(X0-opts.Phi,'fro')/nrmPhi;
   errs_X=err_X;
   eerrs_X=max(max(abs(opts.Phi-X0)./opts.Phi));
end

itn=1;
Nyx=Y0*X0; Nxy=X0*Y0;

if opts.meth==1
   v1wk=max(u1-Nyx*u1,0);
   v2wk=max(u2-Nxy*u2,0);
else % opts.meth==2
   tmp1=E0*u1; tmp2=F0*u2;
   v1wk=w0(1:m)+tmp1+Y0*(tmp2+w0(m+1:mpn));
   v2wk=w0(m+1:mpn)+tmp2+X0*(tmp1+w0(1:m));
end

for i=1:m
    Nyx(i,i)=0.0;
end
for i=1:n
    Nxy(i,i)=0.0;
end

EYX=geMx(Nyx,u1,v1wk,E0,2);
FXY=geMx(Nxy,u2,v2wk,F0,2);
E1=EYX*E0;
F1=FXY*F0;

Xdiff1=FXY*X0*E0; Ydiff1=EYX*Y0*F0;
Y1=Y0+Ydiff1; X1=X0+Xdiff1; 
abserr_upd1=max(max(max(Xdiff1)),max(max(Ydiff1)));
err_upd1=max(max(max(Xdiff1./X1)),max(max(Ydiff1./Y1)));
errs_upd=[abserr_upd0 1; abserr_upd1 err_upd1];

if opts.meth==2
   w1(1:m)=w0(1:m)+EYX*(w0(1:m)+Y0*w0(m+1:mpn));
   w1(m+1:mpn)=w0(m+1:mpn)+FXY*(X0*w0(1:m)+w0(m+1:mpn));
end

if opts.res==1, 
   for i=1:n
       RR(i,:)=diagW(m+i)*X1(i,:);
   end
   for j=1:m
       RR(:,j)=RR(:,j)+X1(:,j)*diagW(j);
   end
   RL=X1*D*X1+N_A*X1+X1*N_B+C;
   
   ResX=RL-RR;
   eres=max(max(abs(ResX)./RR));


   nrmX=norm(X1,'fro'); nrmXp=norm(X1,1);
   res=norm(ResX,'fro')/(nrmX*(nrmXp*nrmD+nrmA+nrmB)+nrmC);
   
   ress=[ress; res eres];
end 

if opts.exact ==1,
   err_X=norm(X1-opts.Phi,'fro')/nrmPhi;
   errs_X=[errs_X; err_X];
   eerrs_X=[eerrs_X; max(max(abs(opts.Phi-X1)./opts.Phi))];
end

tmpX=Xdiff0-Xdiff1;
tmpY=Ydiff0-Ydiff1;
err_stop=max(max(max((Xdiff1.*Xdiff1)./(tmpX.*X1+1e-250))),max(max((Ydiff1.*Ydiff1)./(tmpY.*Y1+1e-250))));
errs_stop=err_stop;
tol=1;
while tol>=opts.tol && itn<opts.maxitn
%while ( (err_stop > opts.tol) || (abserr_upd0 < abserr_upd1) ) && itn<opts.maxitn

    if (abserr_upd0 < abserr_upd1) && (itn>4)
        % monotonicity in update is lost after 4 iterations; possibily convergence but have to check residual to make sure  
       for i=1:n
           RR(i,:)=diagW(m+i)*X1(i,:);
       end
       for j=1:m
           RR(:,j)=RR(:,j)+X1(:,j)*diagW(j);
       end
       RL=X1*D*X1+N_A*X1+X1*N_B+C;
       
       ResX=RL-RR;
       ereswk=max(max(abs(ResX)./RR));
       if ereswk <= sqrt(eps)
          break
       end
    end
    
    E0=E1; F0=F1; X0=X1; Y0=Y1; w0=w1;
    Xdiff0=Xdiff1; Ydiff0=Ydiff1;
    
    Nyx=Y0*X0; Nxy=X0*Y0;
    
    if opts.meth==1
       v1wk=max(u1-Nyx*u1,0);
       v2wk=max(u2-Nxy*u2,0);
    else % opts.meth==2
       tmp1=E0*u1; tmp2=F0*u2;
       v1wk=w0(1:m)+tmp1+Y0*(tmp2+w0(m+1:mpn));
       v2wk=w0(m+1:mpn)+tmp2+X0*(tmp1+w0(1:m));
    end
    
    for i=1:m
        Nyx(i,i)=0.0;
    end
    for i=1:n
        Nxy(i,i)=0.0;
    end
    
    EYX=geMx(Nyx,u1,v1wk,E0,2);
    FXY=geMx(Nxy,u2,v2wk,F0,2);
    E1=EYX*E0;
    F1=FXY*F0;
    
    Xdiff1=FXY*X0*E0; Ydiff1=EYX*Y0*F0;
    Y1=Y0+Ydiff1; X1=X0+Xdiff1; 
    ERX=Xdiff1./(X0+1e-300*ones(n,m));
   ERY=Ydiff1./(Y0+1e-300*ones(m,n));
   tol=max(max(max(abs(ERX))),max(max(abs(ERY))));
    
    abserr_upd0=abserr_upd1;
    abserr_upd1=max(max(max(Xdiff1)),max(max(Ydiff1)));  
    err_upd1=max(max(max(Xdiff1./X1)),max(max(Ydiff1./Y1)));
    errs_upd=[errs_upd; abserr_upd1 err_upd1];  
    
    if opts.meth==2
       w1(1:m)=w0(1:m)+EYX*(w0(1:m)+Y0*w0(m+1:mpn));
       w1(m+1:mpn)=w0(m+1:mpn)+FXY*(X0*w0(1:m)+w0(m+1:mpn));
    end
    
    if opts.res==1, 
       for i=1:n
           RR(i,:)=diagW(m+i)*X1(i,:);
       end
       for j=1:m
           RR(:,j)=RR(:,j)+X1(:,j)*diagW(j);
       end
       RL=X1*D*X1+N_A*X1+X1*N_B+C;
       
       ResX=RL-RR;
       eres=max(max(abs(ResX)./RR));
    
    
       nrmX=norm(X1,'fro'); nrmXp=norm(X1,1);
       res=norm(ResX,'fro')/(nrmX*(nrmXp*nrmD+nrmA+nrmB)+nrmC);
       
       ress=[ress; res eres];
    end 
    
    if opts.exact == 1,
       err_X=norm(X1-opts.Phi,'fro')/nrmPhi;
       errs_X=[errs_X; err_X];
       eerrs_X=[eerrs_X; max(max(abs(opts.Phi-X1)./opts.Phi))];
    end
    
    tmpX=Xdiff0-Xdiff1;
    tmpY=Ydiff0-Ydiff1;
    err_stop=max(max(max((Xdiff1.*Xdiff1)./(tmpX.*X1+1e-300))), ...
        max(max((Ydiff1.*Ydiff1)./(tmpY.*Y1+1e-300))));
    errs_stop=[errs_stop; err_stop];
    
    itn=itn+1; 
    
end

info.itn=itn;
info.upd=errs_upd;
info.stop=errs_stop;
if itn>=opts.maxitn
   info.succ=0;
end

if opts.res==1
   info.res=ress;
end

if opts.exact==1
   info.err_X=errs_X;
   info.eerr_X=eerrs_X;
end

% for use by accADD_L
if opts.meth==2
   info.tw=F1*u2+w1(m+1:m+n);
end

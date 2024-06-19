function [accP,ZZ,Lambda,k]=newJSCaccADDA(AI,BI,CI,DI,n,m,nb,u,v,x1,x2,x3,x4)

p=nb;
N=m+n-p;
A=BI';B=AI';
C=CI';D=DI';
diagA11=diag(A(1:nb,1:nb));
diagB11=diag(B(1:nb,1:nb));
diagC11=diag(C(1:nb,1:nb));
diagD11=diag(D(1:nb,1:nb));

II0=find(diagD11(1:nb)==0); 
II1=find(diagD11(1:nb)>0);

Dt=(diagA11(II1)+diagB11(II1)).^2-4*diagC11(II1).*diagD11(II1);

T0=-diagC11(II0)./(diagA11(II0)+diagB11(II0));
T1=max( -diagB11(II1)./diagD11(II1),  ( (-(diagA11(II1)+diagB11(II1)))+sqrt(Dt) )./(2*diagD11(II1)) );

Lambda0=zeros(p,1);
Lambda0(II0)=T0; Lambda0(II1)=T1;

Lambda=1.01*Lambda0;

AO=A; AO(1:p,1:p)=AO(1:p,1:p)+diag(Lambda)*D(1:p,1:p);
BO=B; BO(1:p,1:p)=BO(1:p,1:p)+D(1:p,1:p)*diag(Lambda);
CO=C; CO(1:p,1:p)=CO(1:p,1:p)+diag(Lambda)*D(1:p,1:p)*diag(Lambda)+A(1:p,1:p)*diag(Lambda)+diag(Lambda)*B(1:p,1:p);

%--------- only works for II0=\emptyset

u0=((diagA11+Lambda.*diagD11)' .* u(m-p+1:m) )./ (diagD11');
v0=v(m-p+1:m)./(diagD11');
v0_hat=zeros(1,p);

%-----------------Left triplet representation for shifted W

WO=[BO -D; -CO AO];
W.N=diag(diag(WO))-WO;
W.u=[u0 u];
W.v=[v0 v(1:m-p) v0_hat v(m+1:N)];

opts.tol=x1;  %  tolerance for stopping DA iteration
opts.maxitn=x2;  % max number of iterations allowed
opts.res=x3;
opts.exact=x4;
%opts.Phi=double(Xacc);


opts.nrms=[norm(AO,1) norm(BO,1) norm(CO,'fro') norm(D,1)];
opts.meth=2;

[X1, Y1, info1]=accADDA_Ldouble(m, n, W, opts);
k=info1.itn;
ZZ=Y1';
% left triplet representation for Z=A-Phi*D
Z.u=W.u(m+1:m+n);
Xnnsft1=X1; Xnnsft1(1:p,1:p)=Xnnsft1(1:p,1:p)-diag(Lambda);
Z.N=A-Xnnsft1*D; Z.N=diag(diag(Z.N))-Z.N;
Z.v=[D(1:p,1:p)*info1.tw(1:p); (D(p+1:m,p+1:n))'*info1.tw(p+1:m)+v(m+1:N)'];
Z.v1=Z.v';
Z1=diag((Z.v1+Z.u*Z.N)./Z.u)-Z.N;  % A-Phi*D
accG=-Z1';
accGamma=X1(1:n,p+1:m)';
accP=[accG(1:p,1:n);accGamma];







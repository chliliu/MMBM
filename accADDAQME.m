function [Gammaplus,Gplus,Gammaminus,Gminus,info]=accADDAQME(nb,e,F,opts)
p=nb;
f=(F.N*F.u+F.v)./F.u;
TF=diag(f)-F.N;
alpha0=max((e+sqrt(e.*e+4*f'))/2);
beta0=max((-e+sqrt(e.*e+4*f'))/2);
E=diag(e);
eta=0.9;
alpha=eta*(1/alpha0);
beta=eta*(1/beta0);
ab=1/(alpha+beta);
E0=ab*(eye(p)-alpha*E-alpha^2*TF);
F0=ab*(eye(p)+beta*E-beta^2*TF);
X0N=ab*alpha*beta*F.N; %offdiag(X0)
Y0N=X0N;
XY0N=X0N+Y0N;% offdiag(X0+Y0)

w1=alpha*F.v;
w2=beta*F.v;
X0v=F0*F.u+w2; %-X0*F.u=X0v
Y0v=E0*F.u+w1;% -Y0*F.u=Y0v

XY0v=X0v+Y0v;% (-X0-Y0)*F.u=XY0v


X0=X0N-diag((X0v+X0N*F.u)./F.u);
Y0=Y0N-diag((Y0v+Y0N*F.u)./F.u);
nI=norm(eye(p),"fro");
nE=norm(E,"fro");
nF=norm(TF,"fro");
tol=1;
itn=0;
r=1;
Xone=1e-300*ones(p,p)+(-2*1e-300)*eye(p);
info.errX=zeros(1,opts.maxitn);
info.EerrX=zeros(1,opts.maxitn);
info.errY=zeros(1,opts.maxitn);
info.EerrY=zeros(1,opts.maxitn);
while ((tol>opts.tol)||(r==0)) && itn<opts.maxitn
   itn=itn+1;
   B=geMx(XY0N,F.u,XY0v,[E0;F0],2);%solve B(-X0-Y0)=[E0;F0], 
   %B=GTH(XY0N,F.u,XY0v,[E0;F0]);
   EG=B(1:p,1:p);% EG=E0(-X0-Y0)^{-1}
   FG=B(p+1:2*p,1:p);%FG=F0(-X0-Y0)^{-1}
   E1=EG*E0;
   F1=FG*F0;
   w=w1+w2;
   w1=w1+EG*w;
   w2=w2+FG*w;

   ErY=EG*F0;
   ErX=FG*E0;
   X1N=X0N+(ErX-diag(diag(ErX)));% offdiag(X1)
   Y1N=Y0N+(ErY-diag(diag(ErY)));% offdiag(Y1)
   X1v=F1*F.u+w2;%-X1*u=X1v
   Y1v=E1*F.u+w1;% -Y1*u=Y1v
   X1=X1N-diag((X1v+X1N*F.u)./F.u);
   Y1=Y1N-diag((Y1v+Y1N*F.u)./F.u);
   ERX=ErX./(X0+Xone);
   ERY=ErY./(Y0+Xone);
   tol=max(max(max(abs(ERX))),max(max(abs(ERY))));
   if opts.exact==1
      info.errX(itn)=norm(X1-opts.Phi,1)/norm(opts.Phi,1);
      info.errY(itn)=norm(Y1-opts.Psi,1)/norm(opts.Psi,1);
      info.EerrX(itn)=max(max(abs((X1-opts.Phi)./(opts.Phi+Xone))));
      info.EerrY(itn)=max(max(abs((Y1-opts.Psi)./(opts.Psi+Xone))));
   end
   E0=E1; F0=F1;   X0N=X1N; Y0N=Y1N; X0=X1; Y0=Y1;
   XY0N=X0N+Y0N;
   XY0v=X1v+Y1v;
   if tol<=opts.tol
       nX0=norm(X0,"fro");
       nY0=norm(Y0,"fro");
       NResX=norm(X0*X0+E*X0-TF,"fro")/(nX0^2*nI+nX0*nE+nF);
       NResY=norm(Y0*Y0-E*Y0-TF,"fro")/(nY0^2*nI+nY0*nE+nF);
       if NResX<opts.tol && NResY<opts.tol
           r=1;
       else 
           r=0;
       end
   end
end
info.itn=itn;
if itn>=opts.maxitn
   info.succ=0;
else 
    info.succ=1;
end

Gammaminus=[];
Gammaplus=[];
Gminus=X0;
Gplus=Y0;










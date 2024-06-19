function [Gammaplus,Gplus,Gammaminus,Gminus,info]=ADDAQME(nb,e,F,opts)
p=nb;
f=diag(F);
alpha0=max((e+sqrt(e.*e+4*f'))/2);
beta0=max((-e+sqrt(e.*e+4*f'))/2);
E=diag(e);
eta=0.9;
alpha=eta*(1/alpha0);
beta=eta*(1/beta0);

TMP=[beta*eye(p) -alpha*eye(p);eye(p)+beta*E eye(p)-alpha*E]\[eye(p) -eye(p);-alpha*F -beta*F];
E0=TMP(1:p,1:p); Y0=TMP(1:p,p+1:2*p);
X0=TMP(p+1:2*p,1:p);F0=TMP(p+1:2*p,p+1:2*p);
 nI=norm(eye(p),"fro");
 nE=norm(E,"fro");
 nF=norm(F,"fro");
tol=1;
itn=0;
r=1;
Xone=1e-300*ones(p,p)+(-2*1e-300)*eye(p);
info.errX=zeros(1,opts.maxitn);
info.EerrX=zeros(1,opts.maxitn);
info.errY=zeros(1,opts.maxitn);
info.EerrY=zeros(1,opts.maxitn);
%while ((tol>opts.tol)||(r==0)) && itn<opts.maxitn
while itn<opts.accitn
    itn=itn+1;
    EXY=E0/(-X0-Y0);
    FXY=F0/(-X0-Y0);
    ErX=FXY*E0;
    ErY=EXY*F0;
    E1=EXY*E0;
    F1=FXY*F0;
    X1=X0+ErX;
    Y1=Y0+ErY;
%    ERX=ErX./(X0+Xone);
%    ERY=ErY./(Y0+Xone);
     tol=max(norm(ErX,1)/norm(X0,1),norm(ErY,1)/norm(Y0,1));
    E0=E1;
    F0=F1;
    X0=X1;
    Y0=Y1;
  if tol<=opts.tol
       nX0=norm(X0,"fro");
       nY0=norm(Y0,"fro");
       NResX=norm(X0*X0+E*X0-F,"fro")/(nX0^2*nI+nX0*nE+nF);
       NResY=norm(Y0*Y0-E*Y0-F,"fro")/(nY0^2*nI+nY0*nE+nF);
       if NResX<opts.tol && NResY<opts.tol
           r=1;
       else 
           r=0;
       end
  end
   if opts.exact==1
       info.errX(itn)=norm(X1-opts.Phi,1)/norm(opts.Phi,1);
       info.errY(itn)=norm(Y1-opts.Psi,1)/norm(opts.Psi,1);
       info.EerrX(itn)=max(max(abs((X1-opts.Phi)./(opts.Phi+Xone))));
       info.EerrY(itn)=max(max(abs((Y1-opts.Psi)./(opts.Psi+Xone))));
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


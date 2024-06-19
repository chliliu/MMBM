function [Gammaplus,Gplus,Gammaminus,Gminus,info]=ADDAMARE(nu,nd,F,opts)
m=nu;
n=nd;
N=m+n;
Fuu=F(1:nu,1:nu);Fud=F(1:nu,nu+1:N);
Fdu=F(nu+1:N,1:nu);Fdd=F(nu+1:N,nu+1:N);
alpha0=max(diag(Fdd));
beta0=max(diag(Fuu));
eta=0.9;
alpha=eta*(1/alpha0);
beta=eta*(1/beta0);
TMP=[beta*Fdd+eye(n), alpha*Fdu; beta*Fud, alpha*Fuu+eye(m)] \ [eye(n)-alpha*Fdd, -beta*Fdu; -alpha*Fud, eye(m)-beta*Fuu];

E0=TMP(1:n,1:n);     Y0=TMP(1:n,n+1:m+n);
X0=TMP(n+1:m+n,1:n); F0=TMP(n+1:m+n,n+1:m+n);
tol=1;
itn=0;
AX=Fuu;BX=Fdd;CX=-Fud;DX=-Fdu;
AY=Fdd;BY=Fuu;CY=-Fdu;DY=-Fud;
nAX=norm(AX,"fro");nBX=norm(BX,"fro");nCX=norm(CX,"fro");nDX=norm(DX,"fro");
nAY=norm(AY,"fro");nBY=norm(BY,"fro");nCY=norm(CY,"fro");nDY=norm(DY,"fro");
r=1;
info.errX=zeros(1,opts.maxitn);
info.EerrX=zeros(1,opts.maxitn);
info.errY=zeros(1,opts.maxitn);
info.EerrY=zeros(1,opts.maxitn);
%while ((tol>opts.tol)||(r==0)) && itn<opts.maxitn
while itn<opts.accitn
   itn=itn+1;
   EYX=E0/(eye(n)-Y0*X0);
   FXY=F0/(eye(m)-X0*Y0);
   Xdiff=FXY*X0*E0;
   Ydiff=EYX*Y0*F0;
   E1=EYX*E0;
   F1=FXY*F0;
   X1=X0+Xdiff;
   Y1=Y0+Ydiff;
%    ERX=Xdiff./(X0+1e-300*ones(m,n));
%    ERY=Ydiff./(Y0+1e-300*ones(n,m));
   tol=max(norm(ErX,1)/norm(X0,1),norm(ErY,1)/norm(Y0,1));
   E0=E1;
   F0=F1;
   X0=X1;
   Y0=Y1;
   if tol<=opts.tol
       nX0=norm(X0,"fro");
       nY0=norm(Y0,"fro");
       NResX=norm(X0*DX*X0-X0*BX-AX*X0+CX,"fro")/(nX0^2*nDX+nX0*nBX+nX0*nAX+nCX);
       NResY=norm(Y0*DY*Y0-Y0*BY-AY*Y0+CY,"fro")/(nY0^2*nDY+nY0*nBY+nY0*nAY+nCY);
       if NResX<opts.tol && NResY<opts.tol
           r=1;
       else 
           r=0;
       end
   end
    if opts.exact==1
       info.errX(itn)=norm(X1-opts.Phi,1)/norm(opts.Phi,1);
       info.errY(itn)=norm(Y1-opts.Psi,1)/norm(opts.Psi,1);
       info.EerrX(itn)=max(max(abs((X1-opts.Phi)./(opts.Phi+1e-300*ones(m,n)))));
       info.EerrY(itn)=max(max(abs((Y1-opts.Psi)./(opts.Psi+1e-300*ones(n,m)))));
    end 
end
info.itn=itn;
if itn>=opts.maxitn
   info.succ=0;
else 
    info.succ=1;
end

Gammaplus=Y1;
Gplus=-Fuu-Fud*Y1;
Gammaminus=X1;
Gminus=-Fdd-Fdu*X1;
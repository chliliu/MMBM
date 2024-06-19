function [Gammaplus,Gplus,Gammaminus,Gminus,info]=accADDAMARE(nu,nd,F,opts)
m=nu;
n=nd;
N=m+n;
diagF=(F.N*F.u+F.v)./F.u;
Fuu=-F.N(1:nu,1:nu)+diag(diagF(1:nu));Fud=-F.N(1:nu,nu+1:N);
Fdu=-F.N(nu+1:N,1:nu);Fdd=-F.N(nu+1:N,nu+1:N)+diag(diagF(nu+1:N));
alpha0=max(diag(Fdd));
beta0=max(diag(Fuu));
eta=0.9;
alpha=eta*(1/alpha0);
beta=eta*(1/beta0);

vu=F.v(1:nu);
vd=F.v(nu+1:N);
uu=F.u(1:nu);
ud=F.u(nu+1:N);

G22v=[1/beta*ud;1/alpha*uu]+[vd;vu];
G22=[eye(n)+beta*Fdd alpha*Fdu;beta*Fud eye(m)+alpha*Fuu];
G22N=diag(diag(G22))-G22;
G22u=[1/beta*ud;1/alpha*uu];
H22=[alpha*Fdd-eye(n) beta*Fdu;alpha*Fud beta*Fuu-eye(m)];

T=geMx(G22N,G22u,G22v,-H22,1);
E0=T(1:nd,1:nd);F0=T(nd+1:N,nd+1:N);
X0=T(nd+1:N,1:nd);Y0=T(1:nd,nd+1:N);

w=geMx(G22N,G22u,G22v,[vd;vu],1);%G22^{-1}[vd;vu]
w10=(alpha+beta)*w(1:n);
w20=(alpha+beta)*w(n+1:N);

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
while ((tol>opts.tol)||(r==0)) && itn<opts.maxitn
   itn=itn+1;
   G10=eye(nd)-Y0*X0;
   G20=eye(nu)-X0*Y0;
   G10N=diag(diag(G10))-G10;
   G20N=diag(diag(G20))-G20;

   u40=w20+F0*uu+X0*(w10+E0*ud);%G2\one=u40
   u30=w10+E0*ud+Y0*(w20+F0*uu);%G1\one=u30
   EG1=geMx(G10N,ud,u30,E0,2);
   FG2=geMx(G20N,uu,u40,F0,2);
   E1=EG1*E0;
   F1=FG2*F0;
   FG2X=FG2*X0;
   EG1Y=EG1*Y0;
   ErX=FG2X*E0;
   ErY=EG1Y*F0;
   X1=X0+ErX;
   Y1=Y0+ErY;
   ERX=ErX./(X0+1e-300*ones(m,n));
   ERY=ErY./(Y0+1e-300*ones(n,m));
   tol=max(max(max(abs(ERX))),max(max(abs(ERY))));
   w1=w10+EG1*(w10+Y0*w20);
   w2=w20+FG2*(X0*w10+w20); 
   E0=E1;
   F0=F1;
   X0=X1;
   Y0=Y1;
   w10=w1;
   w20=w2;
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
Gplusv=vu+(-Fud)*(w10+E0*ud);
GplusN=Gplus-diag(diag(Gplus));
Gplus=GplusN-diag((GplusN*uu+Gplusv)./uu);
Gammaminus=X1;
Gminus=-Fdd-Fdu*X1;
Gminusv=vd+(-Fdu)*(w20+F0*uu);
GminusN=Gminus-diag(diag(Gminus));
Gminus=GminusN-diag((GminusN*ud+Gminusv)./ud);

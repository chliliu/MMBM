function [Gammaplus,Gplus,Gammaminus,Gminus,info]=ADDAMMBM(nu,nb,nd,e,F,opts)
%
% Alternating Directional Doubling Algorithm for LQME
%
% W[I_m\\ Gammaplus](Gplus)^2-K[I_m\\ Gammaplus]Gplus-F[I_m\\Gammaplus]=0......(eq-1)
% 
% W[Gammaminus\\I_n](Gminus)^2+K[Gammaminus\\I_n]Gminus-F[Gammaminus\\I_n]=0...(eq-2)
%
% where 
%             nu  nb   nd            nu   nb   nd           nu  nb   nd
%           [ O          ]  nu     [I_nu         ] nu      [Fuu Fub Fud]nu
%         W=[    I_nb    ]  nb , K=[      E      ] nb  , F=[Fbu Fbb Fbd]nb
%           [           O]  nd     [        -I_nd] nd      [Fdu Fdb Fdd]nd
% E is a diognal matrix and F is a nonsingular or irredicuble singular M-matrix
%                nb           nd                         nu           nb
%Gammaminus=[Gammaminus1   Gammaminus2] nu,  Gammaplus=[Gammaplus1  Gammaplus2] nd
%         nb         nd                    nu         nb
%Gminus=[Gminus11  Gminus12] nb     Gplus=[Gplus11  Gplus12]  nu
%       [Gminus21  Gminus22] nd           [Gplus21  Gplus22]  nb
%---------------------------------------------------
%Using the regrouping technique, (eq-1) can be transformed into an ARE
%              XDX-XB-AX+C=0.......(eq-1')    
%where A=[E  Fbd]  B=[       ]    C=[-Fbb -Fbu]   D=[eye(nb)     ]
%        [   Fdd]    [Fub Fuu]      [-Fdb -Fdu]     [        -Fud] 
% (eq-2) can be transformed into an ARE
%              XDX-XB-AX+C=0.......(eq-2')    
%where A=[-E Fbu]  B=[       ]    C=[-Fbb -Fbd]   D=[eye(nb)     ]
%        [   Fuu]    [Fdb Fdd]      [-Fub -Fud]     [        -Fdu]
%
% Usage:
%
%   [Gammaplus,Gplus,Gammaminus,Gminus,info]=ADDAMMBM(nu,nb,nd,e,F,opts) 
%
%---------------------------------------------------
%
% Input
%
%       nu     positive integer to partition F
%
%       nb     positive integer to partition F
%
%       nd     positive integer to partition F
%       
%       e      a row vector with nb entries and E=diag(e)
%       F      a nonsingular or an irreducible singular M-matrix
%       opts  (structure) 
%             opts.tol        relative entrywise tolerence
%             opts.maxitn:    number of maximal iteration steps, 
%                             default 40;

%             opts.exact      = 1  exact Phi and Psi are given
%                             = 0  no Phi and Psi are provided
%opts.Phi=[Gminus11       Gminus12]
%         [Gammaminus1 Gammaminus2]  the exact solution to (eq-2'), provided only when opts.exact = 1.
%                              
%opts.Psi=[Gplus22        GPlus21]
%         [Gammaplus2  Gammaplus1]   the exact solution to (eq-1'), provided only when opts.exact = 1.
%
% Output
%
%       Gammaplus    (nd by nu+nb) sub-probability matrix(or probability matrix)
%       Gplus        (nu+nb by nu+nb) generator matrix
%       Gammaminus   (nu by nd+nb) sub-probability matrix(or probability matrix)
%       Gminus       (nd+nb by nd+nb) generator matrix
%
%
%       info  (optional) struct
%             info.itn  number of iterations taken
%                        
%             info.succ 1  successful exit after satisfying the stopping tolerance
%                       0  exit after the maxitn is reached.
%             info.errX only referenced when opts.exact=1 for recording relative errors between eact approximation X and Phi
%                             norm(X-Phi,1)
%                         ---------------------
%                              norm(Phi,1)
%             info.errY only referenced when opts.exact=1 for recording relative errors between eact approximation Y and Psi
%                             norm(Y-Psi,1)
%                         ---------------------
%                              norm(Psi,1)
%            info.EerrX only referenced when opts.exact=1 for recording entrywise relative errors between eact approximation X and Phi
%                             (X-Phi)_{i,j}
%                 max_{i,j} ---------------------
%                              (Phi)_{i,j}
%            info.EerrY only referenced when opts.exact=1 for recording entrywise relative errors between eact approximation Y and Psi
%                             (Y-Psi)_{i,j}
%                 max_{i,j} ---------------------
%                              (Psi)_{i,j}
%----------------------------------------------------
m=nb+nu;
n=nb+nd;
N=nb+nu+nd; 
%%%%%%%%%% S=S_b%%%%%%%%%%%%%%%%%%%%%%
if (nu==0) && (nd==0) && (nb~=0)
       [Gammaplus,Gplus,Gammaminus,Gminus,info]=ADDAQME(nb,e,F,opts);
       return
end
%%%%%%%%%% S=S_d\cup S_u%%%%%%%%%%%%%%%%%%%%%%%%%
% (eq-1) and (eq-2) are reduced to MAREs
if(nu~=0) && (nd~=0) && (nb==0) 
  [Gammaplus,Gplus,Gammaminus,Gminus,info]=ADDAMARE(nu,nd,F,opts);
  return
end
%%%%%%%%%%%S=S_b\cup S_d%%%%%%%%%%%%%%%%%%%%%%%%%
if  (nu==0) && (nd~=0) && (nb~=0) 
        Fuu=[];Fub=[];Fud=[];
        Fbu=[];Fbb=F(1:nb,1:nb);Fbd=F(1:nb,nb+1:n);
        Fdu=[];Fdb=F(nb+1:n,1:nb);Fdd=F(nb+1:n,nb+1:n);
        f=diag(Fbb);
      alpha0=max(max((e+sqrt(e.*e+4*f'))/2),max(diag(Fdd)));
      beta0=max((-e+sqrt(e.*e+4*f'))/2);
end
%%%%%%%%%%%%S=S_u\cup S_b %%%%%%%%%%%%%%%%%%%%%%%%%
if (nu~=0) && (nd==0) && (nb~=0)
        Fuu=F(1:nu,1:nu);Fub=F(1:nu,nu+1:m);Fud=[];
        Fbu=F(nu+1:m,1:nu);Fbb=F(nu+1:m,nu+1:m);Fbd=[];
        Fdu=[];Fdb=[];Fdd=[];
        f=diag(Fbb);
      alpha0=max((e+sqrt(e.*e+4*f'))/2);
      beta0=max(max((-e+sqrt(e.*e+4*f'))/2), max(diag(Fuu)));
end
%%%%%%%%%%%%%S=S_u\cup S_b\cup S_dS%%%%%%%%%%%%%%%%%%%%%%%%%
 if (nb~=0)&&(nu~=0)&&(nd~=0) 
       Fuu=F(1:nu,1:nu);Fub=F(1:nu,nu+1:m);Fud=F(1:nu,m+1:N);
        Fbu=F(nu+1:m,1:nu);Fbb=F(nu+1:m,nu+1:m);Fbd=F(nu+1:m,m+1:N);
        Fdu=F(m+1:N,1:nu);Fdb=F(m+1:N,nu+1:m);Fdd=F(m+1:N,m+1:N);
                f=diag(Fbb);
    alpha0=max(max((e+sqrt(e.*e+4*f'))/2), max(diag(Fdd)));
     beta0=max(max((-e+sqrt(e.*e+4*f'))/2), max(diag(Fuu)));
end
p=nb;
E=diag(e);
eta=0.9;
alpha=eta*(1/alpha0);
beta=eta*(1/beta0);
ab=1/(alpha+beta);
%%%%%%%%%%%%%%%%%%%%%%%initial values%%%%%%%%%%%%%%%%%%%%
G1G12=ab*[alpha*eye(p);beta*eye(p)]*[-beta*Fbd -alpha*Fbu];%G11^{-1}*(-G12)
G22=[eye(n-p)+beta*Fdd alpha*Fdu;beta*Fud eye(m-p)+alpha*Fuu];
H21=[alpha*Fdb beta*Fdb; alpha*Fub beta*Fub];
H22=[alpha*Fdd-eye(n-p) beta*Fdu;alpha*Fud beta*Fuu-eye(m-p)];
W21=G22\(-H21);%G22^{-1}*(-H21)
W22=G22\(-H22);%G22^{-1}*(-H22)
W12=ab*[alpha*eye(p);beta*eye(p)]*[-alpha*Fbd -beta*Fbu]+G1G12*W22;
W11=ab*[eye(p)-alpha*E-alpha^2*Fbb -eye(p)+alpha*E-alpha*beta*Fbb; 
    -eye(p)-beta*E-alpha*beta*Fbb eye(p)+beta*E-beta^2*Fbb]+G1G12*W21;
E0=[W11(1:p,1:p) W12(1:p,1:n-p);W21(1:n-p,1:p) W22(1:n-p,1:n-p)];
F0=[W11(p+1:2*p,p+1:2*p) W12(p+1:2*p,n-p+1:n-p+m-p);W21(n-p+1:n-p+m-p,p+1:2*p) W22(n-p+1:n-p+m-p,n-p+1:n-p+m-p)];
X0=[W11(p+1:2*p,1:p) W12(p+1:2*p,1:n-p);W21(n-p+1:n-p+m-p,1:p) W22(n-p+1:n-p+m-p,1:n-p)];
Y0=[W11(1:p,p+1:2*p) W12(1:p,n-p+1:n-p+m-p);W21(1:n-p,p+1:2*p) W22(1:n-p,n-p+1:n-p+m-p)];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tol=1;
itn=0;
AX=[-E Fbu;zeros(nu,p) Fuu];BX=[zeros(p,p) zeros(p,nd);Fdb Fdd];
CX=-[Fbb Fbd;Fub Fud];DX=[eye(p) zeros(p,nu);zeros(nd,p) -Fdu];
AY=[E Fbd;zeros(nd,p) Fdd];BY=[zeros(p,p) zeros(p,nu); Fub Fuu];
CY=-[Fbb Fbu;Fdb Fdu];DY=[eye(p) zeros(p,nd);zeros(nu,p) -Fud];
nAX=norm(AX,"fro");nBX=norm(BX,"fro");nCX=norm(CX,"fro");nDX=norm(DX,"fro");
nAY=norm(AY,"fro");nBY=norm(BY,"fro");nCY=norm(CY,"fro");nDY=norm(DY,"fro");
r=1;
Xone=1e-300*ones(m,n)+[-2*1e-300*eye(p) zeros(p,n-p);zeros(m-p,p) zeros(m-p,n-p)];
Yone=1e-300*ones(n,m)+[-2*1e-300*eye(p) zeros(p,m-p);zeros(n-p,p) zeros(n-p,m-p)];
info.errX=zeros(1,opts.maxitn);
info.EerrX=zeros(1,opts.maxitn);
info.errY=zeros(1,opts.maxitn);
info.EerrY=zeros(1,opts.maxitn);
%while ((tol>opts.tol)||(r==0)) && itn<opts.maxitn
while itn<opts.accitn
   itn=itn+1;
   Z30=[eye(p) X0(1:p,p+1:n);zeros(m-p,p) X0(p+1:m,p+1:n)];
   Z40=[X0(1:p,1:p) zeros(p,m-p); X0(p+1:m,1:p) -eye(m-p)];
   Z50=[Y0(1:p,1:p) zeros(p,n-p);Y0(p+1:n,1:p) -eye(n-p)];
   Z60=[eye(p) Y0(1:p,p+1:m);zeros(n-p,p) Y0(p+1:n,p+1:m)];
   G10=-Z40-Z30*Y0;
   G20=-Z50-Z60*X0;

   EG2=E0/G20;
   FG1=F0/G10;
   E1=EG2*E0;
   F1=FG1*F0;
   FG1Z3=FG1*Z30;
   EG2Z6=EG2*Z60;
   ErX=FG1Z3*E0;
   ErY=EG2Z6*F0;
   X1=X0+ErX;
   Y1=Y0+ErY;
%    ERX=ErX./(X0+Xone);
%    ERY=ErY./(Y0+Yone);
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
      info.EerrX(itn)=max(max(abs((X1-opts.Phi)./(opts.Phi+Xone))));
      info.EerrY(itn)=max(max(abs((Y1-opts.Psi)./(opts.Psi+Yone))));
  end  
end

info.itn=itn;
if itn>=opts.maxitn
   info.succ=0;
else 
    info.succ=1;
end

if  (nu==0) && (nd~=0) && (nb~=0) % S=S_b\cup S_d
   Gammaminus=[];
   Gminus=zeros(n,n);
   Gminus(1:p,1:n)=X0;
   Gminus(p+1:n,1:n)=[-Fdb -Fdd];
   Gammaplus=Y0(p+1:n,1:m);
   Gplus=zeros(m,m);
   Gplus(1:p,1:m)=Y0(1:p,1:m);
end

if (nu~=0) && (nd==0) && (nb~=0)% S=S_u\cup S_b
   Gammaminus=X0(p+1:m,1:n);
   Gminus=zeros(n,n);
   Gminus(1:p,1:n)=X0(1:p,1:n);
   Gammaplus=[];
   Gplus=zeros(m,m);
   Gplus(1:p,1:m)=Y0(1:p,1:m);
   Gplus(p+1:m,1:m)=[-Fub -Fuu];
   P2=[zeros(p,m-p) eye(p);eye(m-p) zeros(m-p,p)];
   Gplus=P2'*Gplus*P2;      
end

if (nb~=0)&&(nu~=0)&&(nd~=0) % S=S_u\cup S_b\cup S_d
   Gammaminus=X0(p+1:m,1:n);
   Gminus=zeros(n,n);
   Gminus(1:p,1:n)=X0(1:p,1:n);
   Gminus(p+1:n,1:n)=-Fdu*Gammaminus+[-Fdb -Fdd];
   Gammaplus=Y0(p+1:n,1:m);
   Gplus=zeros(m,m);
   Gplus(1:p,1:m)=Y0(1:p,1:m);
   Gplus(p+1:m,1:m)=-Fud*Gammaplus+[-Fub -Fuu];
   P2=[zeros(p,m-p) eye(p);eye(m-p) zeros(m-p,p)];
   Gammaplus=Gammaplus*P2;
   Gplus=P2'*Gplus*P2;
end




    


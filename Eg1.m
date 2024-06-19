Q=-[4 0 0 0 0 -4;                    
   0  (15+10^(-8)) -5 -5 -5 -10^(-8);          
   0 -5 15 -5 -5 0;
   0 -5  -5  15 -5 0;
   0 -5 -5 -5 15 0;
   -4 -1 0 0 0 5];

ru=[];
rd=-1.001;
rb=[1 1 1 -1.001 -1.001];
vb=ones(1,5);
nb=5;nu=0;nd=1;p=nb;
m=nu+nb; n=nb+nd;
N=nb+nu+nd;
V=diag([1./ru 1./vb -1./rd]);
x=zeros(N,1);
TF=-V*Q;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Fuu=[];Fub=[];Fud=[];Fbu=[];
Fbb=TF(1:p,1:p);
Fbd=TF(1:p,p+1:n);Fdu=[];
Fdb=TF(p+1:n,1:p);
Fdd=TF(p+1:n,p+1:n);
e=rb./vb;
E=diag(e);
%%%%%%compute accPhi%%%%%%%%%%%%%%%%%%
DI=[eye(nb) ;zeros(nd,nb)];
BI=[zeros(nb,nb) zeros(nb,nd); Fdb Fdd];
AI=-E;
CI=-[Fbb Fbd];
u=ones(1,N);
v=x';
digits(100);
aAI=vpa(AI);
aBI=vpa(BI);
aCI=vpa(CI);
aDI=vpa(DI);
au=vpa(u);
av=vpa(v);
[accPhivpa,ItnJcmPhi]=JSCaccADDA(aAI,aBI,CI,aDI,n,m,nb,au,av);
accPhi=double(accPhivpa);
%%%%%%compute accPsi%%%%%%%%%%%%%%%%%%
DT=[eye(nb) zeros(nb,nd);zeros(nu,nb) -Fud];
BT=[zeros(nb,nb) zeros(nb,nu); Fub Fuu];
AT=[E Fbd;zeros(nd,nb) Fdd];
CT=-[Fbb Fbu;Fdb Fdu];

v=v*[zeros(nu,nd) zeros(nu,nb) eye(nu);zeros(nb,nd) eye(nb) zeros(nb,nu);eye(nd) zeros(nd,nb) zeros(nd,nu)];
digits(100)
aDT=vpa(DT);
aBT=vpa(BT);
aCT=vpa(CT);
aAT=vpa(AT);
av=vpa(v);
[accPsivpa,ItnJcmPsi]=JSCaccADDA(aAT,aBT,aCT,aDT,m,n,nb,au,av);
accPsi=double(accPsivpa);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
MPhi=max(max(abs(accPhi)));
MPsi=max(max(abs(accPsi)));
mPhi=min(min(abs(accPhi)));
mPsi=min(min(abs(accPsi)));
Mmentry=[MPhi,mPhi,MPsi,mPsi];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%accurate DA%%%%%%%%%%%%%%%
opts.tol=1e-14;
opts.maxitn=70;
opts.exact=1;
opts.Phi=accPhi;
opts.Psi=accPsi;
F.N=-TF+diag(diag(TF));
F.v=x;
F.u=ones(N,1);
[Gammaplus,Gplus,Gammaminus,Gminus,infoacc]=accADDAMMBM(nu,nb,nd,e,F,opts);
aRErr1=infoacc.errX;
aRErr2=infoacc.errY;
aERErr1=infoacc.EerrX;
aERErr2=infoacc.EerrY;
figure(1)
semilogy(aRErr1,'g-+','LineWidth',1.4,'MarkerSize',8);
hold on
semilogy(aRErr2,'k--s','LineWidth',1.4,'MarkerSize',8);
hold on
semilogy(aERErr1,'r--x','LineWidth',1.4,'MarkerSize',8);
hold on 
semilogy(aERErr2,'b--o','LineWidth',1.4,'MarkerSize',8);
hold on
legend('RErr(Phi)','RErr(Psi)','ERErr(Phi)','ERErr(Psi)');
xlabel('k'); title('RErr and ERErr');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%plain DA%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
opts.accitn=infoacc.itn;
[Gammaplus,Gplus,Gammaminus,Gminus,info]=ADDAMMBM(nu,nb,nd,e,TF,opts);
RErr1=info.errX;
RErr2=info.errY;
ERErr1=info.EerrX;
ERErr2=info.EerrY;
figure(2)
semilogy(RErr1,'g-+','LineWidth',1.4,'MarkerSize',8);
hold on
semilogy(RErr2,'k--s','LineWidth',1.4,'MarkerSize',8);
hold on
semilogy(ERErr1,'r--x','LineWidth',1.4,'MarkerSize',8);
hold on 
semilogy(ERErr2,'b--o','LineWidth',1.4,'MarkerSize',8);
hold on
legend('RErr(Phi)','RErr(Psi)','ERErr(Phi)','ERErr(Psi)');
xlabel('k'); title('RErr and ERErr');
%%%%%%%%%%%%%%%%%%%%schur method%%%%%%%%%%%%%%%%%
H= [BI -DI;CI -AI];
m=nb+nu;
n=nb+nd;
[U,T] = schur(H,'real'); % Computing the Schur form of H
ee = ordeig(T);
[es,is] = sort(real(ee),'descend');
sel = zeros(m+n,1);
sel(is(1:n)) = 1;
[US,TS]= ordschur(U,T,sel); % Sorting the Schur form of H
schurPhi= US(n+1:m+n,1:n)/US(1:n,1:n);
ErrPhi=max(max(abs((accPhi-schurPhi)./accPhi)));
P3=[eye(nb) zeros(nb,nb+nu+nd);zeros(nd,2*nb+nu) eye(nd);
    zeros(nb,nb+nu) -eye(nb) zeros(nb,nd);zeros(nu,nb) eye(nu) zeros(nu,nb+nd)];
U1=P3'*U;
[es1,is1] = sort(real(-ee),'descend');
sel = zeros(m+n,1);
sel(is1(1:m)) = 1;
[US1,TS1]= ordschur(U1,-T,sel); % Sorting the Schur form of H
schurPsi=US1(m+1:m+n,1:m)/US1(1:m,1:m);
ErrPsi=max(max(abs((accPsi-schurPsi)./accPsi)));
%%%%%%%%extraction DA%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
DI=[eye(nb) ;zeros(nd,nb)];
BI=[zeros(nb,nb) zeros(nb,nd); Fdb Fdd];
AI=-E;
CI=-[Fbb Fbd];
x1=10^(-14);
x2=70;%maxitn
x3=1;x4=0;
[exPhi,Y1,lambda,ItnPsi]=newJSCaccADDA(AI,BI,CI,DI,n,m,nb,u,x',x1,x2,x3,x4);
Lambda=diag(lambda);
Y1bb=Y1(1:nb,1:nb); Y1bu=Y1(1:nb,nb+1:nb+nu);
Y1db=Y1(nb+1:nb+nd,1:nb); Y1du=Y1(nb+1:nd+nb,nb+1:nb+nu);
Y11=[Y1bb Y1bu;zeros(nu,nb) eye(nu)];
Y12=[Lambda*Y1bb-eye(nb) Lambda*Y1bu;Y1db Y1du];
shiftY=Y12/Y11;
%er=norm(Y1*(Lambda*Lambda-E*Lambda-TF)*Y1+DI-Y1*(Lambda-E)-Lambda*Y1,1)/(norm(Y1,1)^2*norm(Lambda*Lambda-E*Lambda-TF,1)+norm(DI,1)+norm((Lambda-E),1)*norm(Y1,1)+norm(Lambda,1)*norm(Y1,1))
%shiftY=-inv(Y1)+Lambda;
shiftErr=norm(shiftY-accPsi,1)/norm(accPsi,1);
shiftEErr=max(max(abs((shiftY-accPsi)./accPsi)));
exPhiEERr=max(max(abs((exPhi-accPhi)./accPhi)));
%err=norm(YOmega-Y1,1)/norm(YOmega,1)
Er=[aERErr1(infoacc.itn) aERErr2(infoacc.itn) shiftEErr aRErr2(infoacc.itn) shiftErr];
fprintf('maxPhi = %d, minPhi = %d, maxPsi = %d, minPsi = %d.\n', ...
    MPhi, mPhi, MPsi, mPsi);
fprintf('ERErraccPhi = %d, ERErraccPsi = %d,ERErrPhi = %d, ERErrPsi = %d, ERErrexPhi = %d,ERErrexPsi = %d, ERErrschurPhi = %d, ERErrschurPsi = %d.\n', ...
    aERErr1(infoacc.itn), aERErr2(infoacc.itn),ERErr1(info.itn), ERErr2(info.itn),exPhiEERr, shiftEErr, ErrPhi,ErrPsi);



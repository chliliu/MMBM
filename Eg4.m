n=16;nu=0;nd=0;nb=n;
x=[-1 1 10^(-1) 5 10^(-1) -10^(3) -10^(-1) 2];
%%%%%%case1
a=10^(-1);b=10^(8.5); vb=1./[b*ones(n/2,1);10^(6)*ones(n/2,1)];
R=diag([-1/2*x -x]./a);
Q=diag(ones(n-1,1),1);
Q(n,1)=1;
Q=Q-eye(n);
Q(n/2,n)=10^(-2);
Q(n/2,n/2)=Q(n/2,n/2)-Q(n/2,n);
Q(n/2+1,n)=10^(-6);
Q(n/2+1,n/2+1)=Q(n/2+1,n/2+1)-Q(n/2+1,n);
u=ones(n,1);
v=zeros(n,1);
N=-Q+diag(diag(Q));
[LF,UF]=geMLU(N,u,v);
% stationary vector pi of F: pi*F=0.
piQ=zeros(1,n); piQ(n)=1;
for i=n-1:-1:1
    piQ(i)=(-piQ(i+1:n)*LF(i+1:n,i) )/LF(i,i);
end
piQ=piQ/sum(piQ);
mu=piQ*R*u
Sigma=diag(vb);
E=diag(1./vb)*R;
TF=diag(1./vb)*(-Q);
e=diag(E)';
%%%%%%compute accPhi%%%%%%%%%%%%%%%%%%
DI=eye(n);
BI=zeros(n,n);
AI=-E;
CI=-TF;
digits(100);
aAI=vpa(AI);
aBI=vpa(BI);
aCI=vpa(CI);
aDI=vpa(DI);
au=vpa(u');
av=vpa(v');
% x1=10^(-80);
% x2=100;%maxitn
% x3=1;
% x4=0;
[accPhivpa,ItnJcmPhi]=JSCaccADDA(aAI,aBI,CI,aDI,n,n,n,au,av);
%[accPhivpa,newZ,Lambda,ItnJcmPhi]=newJSCaccADDA(aAI,aBI,CI,aDI,n,n,n,au,av,x1,x2,x3,x4);
accPhi=double(accPhivpa);
%YOmega=double(newZ);
%%%%%%compute accPsi%%%%%%%%%%%%%%%%%%
DT=eye(n);
BT=zeros(n,n);
AT=E;
CT=-TF;
digits(100)
aDT=vpa(DT);
aBT=vpa(BT);
aCT=vpa(CT);
aAT=vpa(AT);
au=vpa(u');
av=vpa(v');
[accPsivpa,ItnJcmPsi]=JSCaccADDA(aAT,aBT,aCT,aDT,n,n,n,au,av);
%maxitn=100;
accPsi=double(accPsivpa);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
MPhi=max(max(abs(accPhi)));
MPsi=max(max(abs(accPsi)));
mPhi=min(min(abs(accPhi)));
mPsi=min(min(abs(accPsi)));
Mmentry=[MPhi,mPhi,MPsi,mPsi];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%accurate DA%%%%%%%%%%%%%%%%%%%%%%%%%%%
opts.tol=1e-14;
opts.maxitn=70;
opts.exact=1;
opts.Phi=accPhi;
opts.Psi=accPsi;
F.N=-TF+diag(diag(TF));
F.v=zeros(n,1);
F.u=ones(n,1);
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

% %%%%%%%%%%%%%%%%%%%%%%plain DA%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%extraction DA %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
DI=eye(n);
BI=zeros(n,n);
AI=-E;
CI=-TF;
% opts.tol=1e-14;  %  tolerance for stopping DA iteration
% opts.maxitn=50;  % max number of iterations allowed
% opts.res=1;
% opts.exact=0;
x1=10^(-14);
x2=70;%maxitn
x3=1;x4=0;
[exPhi,Y1,lambda,ItnPsi]=newJSCaccADDA(AI,BI,CI,DI,n,n,n,u',v',x1,x2,x3,x4);
Lambda=diag(lambda);
%er=norm(Y1*(Lambda*Lambda-E*Lambda-TF)*Y1+DI-Y1*(Lambda-E)-Lambda*Y1,1)/(norm(Y1,1)^2*norm(Lambda*Lambda-E*Lambda-TF,1)+norm(DI,1)+norm((Lambda-E),1)*norm(Y1,1)+norm(Lambda,1)*norm(Y1,1))
shiftY=-inv(Y1)+Lambda;
shiftErr=norm(shiftY-accPsi,1)/norm(accPsi,1);
shiftEErr=max(max(abs((shiftY-accPsi)./accPsi)));
exPhiEERr=max(max(abs((exPhi-accPhi)./accPhi)));
%err=norm(YOmega-Y1,1)/norm(YOmega,1)
Er=[aERErr1(infoacc.itn) aERErr2(infoacc.itn) shiftEErr aRErr2(infoacc.itn) shiftErr];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%Schur method%%%%%%%%%%%%%%%%%%%%%%%
H= [BI -DI;CI -AI];
m=n;
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
fprintf('---------------------------------------------------------------------------------------------------------------------\n');
fprintf('when a=10^(-8.5), b=10^(-6) \n');
fprintf('maxPhi = %d, minPhi = %d, maxPsi = %d, minPsi = %d.\n', ...
    MPhi, mPhi, MPsi, mPsi);
fprintf('ERErraccPhi = %d, ERErraccPsi = %d,ERErrPhi = %d, ERErrPsi = %d, ERErrexPhi = %d,ERErrexPsi = %d, ERErrschurPhi = %d, ERErrschurPsi = %d.\n', ...
    aERErr1(infoacc.itn), aERErr2(infoacc.itn),ERErr1(info.itn), ERErr2(info.itn),exPhiEERr, shiftEErr, ErrPhi,ErrPsi);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%case2%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% a=10^(-1); b=10^(-7);vb=1./[b*ones(n/2,1);10^(-3)*ones(n/2,1)];
% R=diag([-1/2*x -x]./a);
% Q=diag(ones(n-1,1),1);
% Q(n,1)=1;
% Q=Q-eye(n);
% Q(n/2,n)=10^(-2);
% Q(n/2,n/2)=Q(n/2,n/2)-Q(n/2,n);
% Q(n/2+1,n)=10^(-6);
% Q(n/2+1,n/2+1)=Q(n/2+1,n/2+1)-Q(n/2+1,n);
% u=ones(n,1);
% v=zeros(n,1);
% N=-Q+diag(diag(Q));
% [LF,UF]=geMLU(N,u,v);
% % stationary vector pi of F: pi*F=0.
% piQ=zeros(1,n); piQ(n)=1;
% for i=n-1:-1:1
%     piQ(i)=(-piQ(i+1:n)*LF(i+1:n,i) )/LF(i,i);
% end
% piQ=piQ/sum(piQ);
% mu=piQ*R*u;
% Sigma=diag(vb);
% E=diag(1./vb)*R;
% TF=diag(1./vb)*(-Q);
% e=diag(E)';
% %%%%%%compute accPhi%%%%%%%%%%%%%%%%%%
% DI=eye(n);
% BI=zeros(n,n);
% AI=-E;
% CI=-TF;
% digits(100);
% aAI=vpa(AI);
% aBI=vpa(BI);
% aCI=vpa(CI);
% aDI=vpa(DI);
% au=vpa(u');
% av=vpa(v');
% [accPhivpa,ItnJcmPhi]=JSCaccADDA(aAI,aBI,CI,aDI,n,n,n,au,av);
% accPhi=double(accPhivpa);
% %YOmega=double(newZ);
% %%%%%%compute accPsi%%%%%%%%%%%%%%%%%%
% DT=eye(n);
% BT=zeros(n,n);
% AT=E;
% CT=-TF;
% digits(100)
% aDT=vpa(DT);
% aBT=vpa(BT);
% aCT=vpa(CT);
% aAT=vpa(AT);
% au=vpa(u');
% av=vpa(v');
% [accPsivpa,ItnJcmPsi]=JSCaccADDA(aAT,aBT,aCT,aDT,n,n,n,au,av);
% %maxitn=100;
% accPsi=double(accPsivpa);
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MPhi=max(max(abs(accPhi)));
% MPsi=max(max(abs(accPsi)));
% mPhi=min(min(abs(accPhi)));
% mPsi=min(min(abs(accPsi)));
% Mmentry=[MPhi,mPhi,MPsi,mPsi];
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%accurate DA%%%%%%%%%%%%%%%%%%%%%%%%%%%
% opts.tol=1e-14;
% opts.maxitn=70;
% opts.exact=1;
% opts.Phi=accPhi;
% opts.Psi=accPsi;
% F.N=-TF+diag(diag(TF));
% F.v=zeros(n,1);
% F.u=ones(n,1);
% [Gammaplus,Gplus,Gammaminus,Gminus,infoacc]=accADDAMMBM(nu,nb,nd,e,F,opts);
% aRErr1=infoacc.errX;
% aRErr2=infoacc.errY;
% aERErr1=infoacc.EerrX;
% aERErr2=infoacc.EerrY;
% figure(3)
% semilogy(aRErr1,'g-+','LineWidth',1.4,'MarkerSize',8);
% hold on
% semilogy(aRErr2,'k--s','LineWidth',1.4,'MarkerSize',8);
% hold on
% semilogy(aERErr1,'r--x','LineWidth',1.4,'MarkerSize',8);
% hold on 
% semilogy(aERErr2,'b--o','LineWidth',1.4,'MarkerSize',8);
% hold on
% legend('RErr(Phi)','RErr(Psi)','ERErr(Phi)','ERErr(Psi)');
% xlabel('k'); title('RErr and ERErr');
% 
% % %%%%%%%%%%%%%%%%%%%%%%%%plain DA%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% opts.accitn=infoacc.itn;
% [Gammaplus,Gplus,Gammaminus,Gminus,info]=ADDAMMBM(nu,nb,nd,e,TF,opts);
% RErr1=info.errX;
% RErr2=info.errY;
% ERErr1=info.EerrX;
% ERErr2=info.EerrY;
% figure(4)
% semilogy(RErr1,'g-+','LineWidth',1.4,'MarkerSize',8);
% hold on
% semilogy(RErr2,'k--s','LineWidth',1.4,'MarkerSize',8);
% hold on
% semilogy(ERErr1,'r--x','LineWidth',1.4,'MarkerSize',8);
% hold on 
% semilogy(ERErr2,'b--o','LineWidth',1.4,'MarkerSize',8);
% hold on
% legend('RErr(Phi)','RErr(Psi)','ERErr(Phi)','ERErr(Psi)');
% xlabel('k'); title('RErr and ERErr');
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % %%%%%%%%%%%%%%%%extraction DA%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DI=eye(n);
% BI=zeros(n,n);
% AI=-E;
% CI=-TF;
% % opts.tol=1e-14;  %  tolerance for stopping DA iteration
% % opts.maxitn=50;  % max number of iterations allowed
% % opts.res=1;
% % opts.exact=0;
% x1=10^(-14);
% x2=70;%maxitn
% x3=1;x4=0;
% [exPhi,Y1,lambda,ItnPsi]=newJSCaccADDA(AI,BI,CI,DI,n,n,n,u',v',x1,x2,x3,x4);
% Lambda=diag(lambda);
% %er=norm(Y1*(Lambda*Lambda-E*Lambda-TF)*Y1+DI-Y1*(Lambda-E)-Lambda*Y1,1)/(norm(Y1,1)^2*norm(Lambda*Lambda-E*Lambda-TF,1)+norm(DI,1)+norm((Lambda-E),1)*norm(Y1,1)+norm(Lambda,1)*norm(Y1,1))
% shiftY=-inv(Y1)+Lambda;
% shiftErr=norm(shiftY-accPsi,1)/norm(accPsi,1);
% shiftEErr=max(max(abs((shiftY-accPsi)./accPsi)));
% exPhiEERr=max(max(abs((exPhi-accPhi)./accPhi)));
% %err=norm(YOmega-Y1,1)/norm(YOmega,1)
% Er=[aERErr1(infoacc.itn) aERErr2(infoacc.itn) shiftEErr aRErr2(infoacc.itn) shiftErr];
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%Schur method%%%%%%%%%%%%%%%%%%%%%%%
% H= [BI -DI;CI -AI];
% m=n;
% [U,T] = schur(H,'real'); % Computing the Schur form of H
% ee = ordeig(T);
% [es,is] = sort(real(ee),'descend');
% sel = zeros(m+n,1);
% sel(is(1:n)) = 1;
% [US,TS]= ordschur(U,T,sel); % Sorting the Schur form of H
% schurPhi= US(n+1:m+n,1:n)/US(1:n,1:n);
% ErrPhi=max(max(abs((accPhi-schurPhi)./accPhi)));
% P3=[eye(nb) zeros(nb,nb+nu+nd);zeros(nd,2*nb+nu) eye(nd);
%     zeros(nb,nb+nu) -eye(nb) zeros(nb,nd);zeros(nu,nb) eye(nu) zeros(nu,nb+nd)];
% U1=P3'*U;
% [es1,is1] = sort(real(-ee),'descend');
% sel = zeros(m+n,1);
% sel(is1(1:m)) = 1;
% [US1,TS1]= ordschur(U1,-T,sel); % Sorting the Schur form of H
% schurPsi=US1(m+1:m+n,1:m)/US1(1:m,1:m);
% ErrPsi=max(max(abs((accPsi-schurPsi)./accPsi)));
% fprintf('---------------------------------------------------------------------------------------------------------------------\n');
% fprintf('when a=10^7, b=10^3\n');
% fprintf('maxPhi = %d, minPhi = %d, maxPsi = %d, minPsi = %d.\n', ...
%     MPhi, mPhi, MPsi, mPsi);
% fprintf('ERErraccPhi = %d, ERErraccPsi = %d,ERErrPhi = %d, ERErrPsi = %d, ERErrexPhi = %d,ERErrexPsi = %d, ERErrschurPhi = %d, ERErrschurPsi = %d.\n', ...
%     aERErr1(infoacc.itn), aERErr2(infoacc.itn),ERErr1(info.itn), ERErr2(info.itn),exPhiEERr, shiftEErr, ErrPhi,ErrPsi);
function [K, gam] = h2AS(GP,nmeas,ncon,sub)
% Distributed H2 control for alpha-heterogneous string interconnected
% systems
% Licensed under GPLv3. See License.txt in the project root for license information.
% Author: Simon Heinke

% Inputs:
% GP: Generalized Plant, with subsystems GP.Sub{i}, a row vector GP.Ord
%       specifying the type of subsystem i, a row vector GP.Nki specifying
%       the number of subsystem type ki and GP.ns specifying the order of
%       the interconnection signals.
% nmeas: row vector containing the number of measured outputs of subsystem ki
% ncon: row vector containing the number of control inputs of subsystem ki
% sub: scalar >=1 specifying suboptimality of the bound (sometimes necessary 
%       due to numerical issues)

% Outputs:
% K: alpha-heterogeneous string interconnected controller
% gam: H2 performance bound

ALPHA=size(GP.Sub,2);
gammasq=sdpvar(1);
ns=GP.ns;

% Define the optimisation variables
for i=1:ALPHA
    Sub=GP.Sub{i};
    if i==1 || i==ALPHA
        nsi=ns;
    else
        nsi=2*ns;
    end
    nt=size(Sub.A,1)-nsi;
    
    YT{i}=sdpvar(nt);
    XT{i}=sdpvar(nt);

    alpha{i}=sdpvar(1);
    beta{i}=sdpvar(1);
end

X11=sdpvar(ns);
X22=sdpvar(ns);
Y11=sdpvar(ns);
Y22=sdpvar(ns);
X12=sdpvar(ns,ns,'full');
Y12=sdpvar(ns,ns,'full');

% Construct the LMIs
LMI=[];
Sum=0;
for i=1:ALPHA
    Sub=GP.Sub{i};
    if i==1
        nsi=ns;
        Mki{i}=[X22 -X12'; -X12 -X11];
        Yki{i}=[Y22 -Y12'; -Y12 -Y11];
    elseif i==ALPHA
        nsi=ns;
        Mki{i}=[X11 X12; X12' -X22];
        Yki{i}=[Y11 Y12; Y12' -Y22];
    else
        nsi=2*ns;
        P11=blkdiag(X11,X22);
        P12=[zeros(ns) X12; -X12' zeros(ns)];
        Z11=blkdiag(Y11,Y22);
        Z12=[zeros(ns) Y12; -Y12' zeros(ns)];
        Mki{i}=[P11 P12; P12' -P11];
        Yki{i}=[Z11 Z12; Z12' -Z11];
    end
    nt=size(Sub.A,1)-nsi;
    nu=ncon(i);
    nd=size(Sub.B,2)-nu;
    ny=nmeas(i);
    nz=size(Sub.C,1)-ny;
    
    Att=Sub.A(1:nt,1:nt);
    Ats=Sub.A(1:nt,nt+1:end);
    Btd=Sub.B(1:nt,1:nd);
    Btu=Sub.B(1:nt,nd+1:end);
    Ast=Sub.A(nt+1:end,1:nt);
    Ass=Sub.A(nt+1:end,nt+1:end);
    Bsd=Sub.B(nt+1:end,1:nd);
    Bsu=Sub.B(nt+1:end,nd+1:end);
    Ctz=Sub.C(1:nz,1:nt);
    Csz=Sub.C(1:nz,nt+1:end);
    Dzd=Sub.D(1:nz,1:nd);
    Dzu=Sub.D(1:nz,nd+1:end);
    Cty=Sub.C(nz+1:end,1:nt);
    Csy=Sub.C(nz+1:end,nt+1:end);
    Dyd=Sub.D(nz+1:end,1:nd);
    
    Nx=null([Cty Csy Dyd]);

    Vx=[     eye(nt)         zeros(nt,nsi)     zeros(nt,nd); ...
              Att               Ats              Btd; ...
          zeros(nsi,nt)       eye(nsi)         zeros(nsi,nd); ... 
                Ast             Ass             Bsd; ...
          zeros(nd,nt)      zeros(nd,nsi)       eye(nd); ... 
               Ctz              Csz             Dzd];

    BlockGamma=blkdiag(-alpha{i}*eye(nd),eye(nz));
    BlockGammaInv=blkdiag(-beta{i}*eye(nd),eye(nz));

    Mx=blkdiag([zeros(nt) XT{i}; XT{i} zeros(nt)],Mki{i},BlockGamma);
    Lx=Nx'*Vx'*Mx*Vx*Nx;

    Ny=null([Btu' Bsu' Dzu']);

    Vy=[      Att'             Ast'            Ctz'; ...
            -eye(nt)       zeros(nt,nsi)     zeros(nt,nz); ...
               Ats'             Ass'             Csz'; ...
           zeros(nsi,nt)       -eye(nsi)       zeros(nsi,nz); ...
               Btd'             Bsd'             Dzd'; ...
            zeros(nz,nt)      zeros(nz,nsi)       -eye(nz)];

    My=blkdiag([zeros(nt) YT{i}; YT{i} zeros(nt)],Yki{i},BlockGammaInv);
    Ly=Ny'*Vy'*My*Vy*Ny;

    Lxy=[XT{i} eye(nt); eye(nt) YT{i}];
    
    LMI=[LMI, Lx<=0, Ly>=0, XT{i}>=0, YT{i}>=0, Lxy>=0, alpha{i}>=0, beta{i}>=0];
    
    Sum=Sum+GP.Nki(i)*trace(Btd'*XT{i}*Btd);
end
  
LMI1=[LMI, Sum<=gammasq];

OPTIONS = sdpsettings('solver', 'sedumi');
OPTIONS.showprogess=1;

% solve the LMIs:
diagnostics = optimize(LMI1,gammasq,OPTIONS)
gammasq=value(gammasq);
if sub~=1
    % solve the LMIs again with a suboptimal bound due to numerical issues
    gammasq=gammasq*sub;
    LMI2=[LMI, Sum<=gammasq];
    diagnostics = optimize(LMI2,[],OPTIONS)
end

% Construct closed loop Lyapunov matrix
for i=1:ALPHA    
    YT{i}=value(YT{i});
    XT{i}=value(XT{i});
    alpha{i}=value(alpha{i});
    beta{i}=value(beta{i});
    nt=size(XT{i},1);
    
    [Qtemp,Rtemp] = qr(eye(nt)-XT{i}*YT{i});
    Xgk=Qtemp;
    Ygk=Rtemp';
    Xk=-Xgk'*YT{i}/Ygk';
    Xk=(Xk+Xk')./2;
    XC{i}=[XT{i} Xgk; Xgk' Xk];
end

% Construct cloosed loop multipliers
Case=1; % Case 1: n_{ij}^K=n_{ij}, Case 2: n_{ij}^K=3n_{ij}.
Mk3=value(Mki{ALPHA});
Yk3=value(Yki{ALPHA});
Mk3(isnan(Mk3))=0;
Yk3(isnan(Yk3))=0;

T=[Mk3 eye(size(Mk3,1)); eye(size(Mk3,2)) Yk3];
NoNeg=sum(eig(T)<0); % Number of negative Eigenvalues
NoPos=size(T,1)-NoNeg;
if NoNeg~=NoPos
    Case=2; % if for any egde NoNeg~=NoPos switch to Case 2
end
   
if Case==1
    [Q,R] = qr(eye(2*ns)-Mk3*Yk3);
    S12=Q;
    R12=R';
    S22=-S12'*Yk3/R12';
    S22=(S22+S22')./2;
    MC=[Mk3 S12; S12' S22];
    TC=[eye(ns) zeros(ns,3*ns); zeros(ns,2*ns) eye(ns) zeros(ns); ...
        zeros(ns) eye(ns) zeros(ns,2*ns); zeros(ns,3*ns) eye(ns)];
    MCT{3}=TC'*MC*TC;
    
    Tp=[zeros(2*ns) eye(2*ns); eye(2*ns) zeros(2*ns)];
    MCT{1}=-Tp*MCT{3}*Tp;
    Tp2=[blkdiag(eye(ns),[zeros(ns) eye(ns)],[zeros(ns,2*ns) eye(ns)],[zeros(ns) eye(ns)]); ...
        zeros(2*ns,4*ns) blkdiag([eye(ns) zeros(ns)],[eye(ns) zeros(ns)]); ...
        blkdiag([zeros(ns) eye(ns)],[zeros(ns) eye(ns)]) zeros(2*ns,4*ns)];
    MCT{2}=Tp2'*blkdiag(MCT{3},-MCT{3})*Tp2;
    
elseif Case==2
    T=Mk3-inv(Yk3);
    [V,D]=eig(T);
    [~,I] = sort(diag(D),'descend');
    V=V(:,I);
    D=D(I,I);
    V=V*sqrt(abs(D));
    D=sign(D);
    p=sum(sum(D>0));
    Vp=V(:,1:p);
    Vm=V(:,p+1:end);
    R1p=sum(sign(eig(Yk3))>0);
    R1m=sum(sign(eig(Yk3))<0);
    p=4*ns-R1p;
    m=4*ns-R1m;
    J=blkdiag(eye(p),-eye(m));
    S22=J;
    S12=[Vp zeros(size(Vp,1),p-size(Vp,2)) Vm zeros(size(Vm,1),m-size(Vm,2))];
    MC=[Mk3 S12; S12' S22];
    nsk=3*ns;
    TC=[eye(ns) zeros(ns,ns+2*nsk); zeros(ns,ns+nsk) eye(ns) zeros(ns,nsk); ...
        zeros(nsk,ns) eye(nsk) zeros(nsk,ns+nsk); zeros(nsk,2*ns+nsk) eye(nsk)];
    MCT{3}=TC'*MC*TC;
    
    Tp=[zeros(ns+nsk) eye(ns+nsk); eye(ns+nsk) zeros(ns+nsk)];
    MCT{1}=-Tp'*MCT{3}*Tp;
    Tp2=[blkdiag(eye(ns),[zeros(nsk,ns) eye(nsk)],[zeros(ns,ns+nsk) eye(ns)],[zeros(nsk) eye(nsk)]); ...
        zeros(ns+nsk,2*ns+2*nsk) blkdiag([eye(ns) zeros(ns)],[eye(nsk) zeros(nsk)]); ...
        blkdiag([zeros(ns) eye(ns)],[zeros(nsk) eye(nsk)]) zeros(ns+nsk,2*ns+2*nsk)];
    MCT{2}=Tp2'*blkdiag(MCT{3},-MCT{3})*Tp2;
 
end
%% Controller Reconstruction

for i=1:ALPHA
    if i==1
        MCTi=MCT{1};
        nsi=ns;
    elseif i==ALPHA
        MCTi=MCT{3};
        nsi=ns;
    else
        MCTi=MCT{2};
        nsi=2*ns;
    end
    if Case==1
        nski=nsi;
    elseif Case==2
        nski=3*nsi;
    end
    
    Sub=GP.Sub{i};
    nt=size(Sub.A,1)-nsi;
    nu=ncon(i);
    nd=size(Sub.B,2)-nu;
    ny=nmeas(i);
    nz=size(Sub.C,1)-ny;
    
    Att=Sub.A(1:nt,1:nt);
    Ats=Sub.A(1:nt,nt+1:end);
    Btd=Sub.B(1:nt,1:nd);
    Btu=Sub.B(1:nt,nd+1:end);
    Ast=Sub.A(nt+1:end,1:nt);
    Ass=Sub.A(nt+1:end,nt+1:end);
    Bsd=Sub.B(nt+1:end,1:nd);
    Bsu=Sub.B(nt+1:end,nd+1:end);
    Ctz=Sub.C(1:nz,1:nt);
    Csz=Sub.C(1:nz,nt+1:end);
    Dzd=Sub.D(1:nz,1:nd);
    Dzu=Sub.D(1:nz,nd+1:end);
    Cty=Sub.C(nz+1:end,1:nt);
    Csy=Sub.C(nz+1:end,nt+1:end);
    Dyd=Sub.D(nz+1:end,1:nd);
    
    MCf11=MCTi(1:nsi+nski,1:nsi+nski);
    MCf12=MCTi(1:nsi+nski,nsi+nski+1:end);
    MCf22=MCTi(nsi+nski+1:end,nsi+nski+1:end);
    
    Btuc=[zeros(nt,nt+nski) Btu; eye(nt) zeros(nt,nski+nu)];
    Bsuc=[zeros(nsi,nt+nski) Bsu; zeros(nski,nt) eye(nski) zeros(nski,nu)];
    Dzuc=[zeros(nz,nt+nski) Dzu];
    V=[Btuc;Bsuc;Dzuc]';

    Ctyc=[zeros(nt) eye(nt); zeros(nski,2*nt); Cty zeros(ny,nt)];
    Csyc=[zeros(nt,nsi+nski); zeros(nski,nsi)  eye(nski); Csy zeros(ny,nski)];
    Dydc=[zeros(nt+nski,nd); Dyd];
    U=[Ctyc Csyc Dydc];

    R11=[Att zeros(nt) Ats zeros(nt,nski); zeros(nt,2*nt+nsi+nski); Ast zeros(nsi,nt) Ass zeros(nsi,nski); zeros(nski,2*nt+nsi+nski)];
    R12=[Btd; zeros(size(Btd)); Bsd; zeros(nski,nd)];
    R21=[Ctz zeros(size(Ctz)) Csz zeros(nz,nski)];
    R22=Dzd;
    R=[R11 R12; R21 R22];

    if alpha{i}*beta{i}>=1
        rho{i}=alpha{i};
    else
        rho{i}=1/beta{i};
    end
    
    nc=size(XC{i},1);
    M11=[zeros(nc,nc+nsi+nski+nd); zeros(nsi+nski,nc) MCf11 zeros(nsi+nski,nd); zeros(nd,nc+nsi+nski) -rho{i}*eye(nd)];
    M12=[XC{i} zeros(nc,nsi+nski+nz); zeros(nsi+nski,nc) MCf12 zeros(nsi+nski,nz); zeros(nd,nc+nsi+nski+nz)];
    M22=[zeros(nc,nc+nsi+nski+nz); zeros(nsi+nski,nc) MCf22 zeros(nsi+nski,nz); zeros(nz,nc+nsi+nski) eye(nz)];
    M=[M11 M12; M12' M22];

    [~, ~, H]=svd(V);
    [~, ~, J]=svd(U);

    N=blkdiag(J',inv(H))*M*blkdiag(J,inv(H'));
    N=(N+N')./2;
    T=H'*R*J;
    
    n1=nt+nski+nu;
    n2=nt+nski+ny;
    T11=T(1:n1,1:n2);
    T12=T(1:n1,n2+1:end);
    T21=T(n1+1:end,1:n2);
    T22=T(n1+1:end,n2+1:end);
    
    n1=size(T21,2);
    n2=size(T12,1);
    n3=size(T12,2);
    n4=size(T21,1);
    Q=[eye(n1) zeros(n1,n2); zeros(n3,n1+n2); zeros(n2,n1) eye(n2); T21 zeros(n4,n2)];

    S=[zeros(n1,n3); eye(n3); T12; T22];

    Ma=Q'*N*Q-Q'*N*S/(S'*N*S)*S'*N*Q;
    Ma=(Ma+Ma')./2;

    n1=nt+nski+nu;
    n2=nt+nski+ny;
    [Vt,Dt] = eig(Ma);
    [~,ind] = sort(diag(Dt));
    Vs = Vt(:,ind);

    Z=Vs(:,1:n2);
    Z1=Z(1:n2,:);
    Z2=Z(n2+1:end,:);
    
    Z11=Z2/Z1;

    V1=V*H;
    V1=V1(:,1:n1);
    U1=U*J;
    U1=U1(:,1:n2);

    Theta=V1'\(Z11-T11)/U1;
    K.Sub{i}.A=Theta(1:nt+nski,1:nt+nski);
    K.Sub{i}.B=Theta(1:nt+nski,nt+nski+1:end);
    K.Sub{i}.C=Theta(nt+nski+1:end,1:nt+nski);
    K.Sub{i}.D=Theta(nt+nski+1:end,nt+nski+1:end);
end
K.Ts=GP.Ts;
K.Ord=GP.Ord;
if Case==1
    K.ns=GP.ns;
elseif Case==2
    K.ns=3*GP.ns;
end
K.Nki=GP.Nki;
gam=sqrt(gammasq);




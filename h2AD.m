function [K, gam] = h2AD(GP,sub)
% Distributed H2 control for alpha-heterogenous decomposable systems
% Licensed under GPLv3. See License.txt in the project root for license information.
% Author: Simon Heinke

% Inputs:
% GP: Generalized Plant, with subsystems GP.Sub{i}, interconnection matrix
%       GP.T, size of the interconnection signals GP.ns and a row vector 
%       GP.Ord specifying the type of subsystem i.
% sub: scalar >=1 specifying suboptimality of the bound (sometimes necessary 
%       due to numerical issues)

% Outputs:
% K: alpha-heterogeneous decomposable controller
% gam: H2 performance bound

ALPHA=size(GP.Sub,2); % # of different subsystems
gammasq=sdpvar(1);
ns=GP.ns;
% DG-scalings
D=sdpvar(ns);
G=sdpvar(ns,ns,'skew');
M1=[D G; G' -D];
Di=sdpvar(ns);
Gi=sdpvar(ns,ns,'skew');
M2=[Di Gi; Gi' -Di];

% Construct the LMIs
LMI=[];
Sum=0;
for i=1:ALPHA
    Sub=GP.Sub{i};
    nt=size(Sub.A,1)-ns;
    nu=GP.Sub{i}.nu;
    nd=size(Sub.B,2)-nu;
    ny=GP.Sub{i}.ny;
    nz=size(Sub.C,1)-ny;
    
    YT{i}=sdpvar(nt);
    XT{i}=sdpvar(nt);
    
    alpha{i}=sdpvar(1);
    beta{i}=sdpvar(1);

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

    Vx=[      eye(nt)         zeros(nt,ns)     zeros(nt,nd); ...
                Att             Ats              Btd; ...
           zeros(ns,nt)       eye(ns)         zeros(ns,nd); ...
                Ast             Ass             Bsd; ...
           zeros(nd,nt)      zeros(nd,ns)       eye(nd); ...
               Ctz              Csz             Dzd]; 
            
    BlockGamma=blkdiag(-alpha{i}*eye(nd),eye(nz));
    BlockGammaInv=blkdiag(-beta{i}*eye(nd),eye(nz));

    Mx=blkdiag([zeros(nt) XT{i}; XT{i} zeros(nt)],M1,BlockGamma);
    Lx=Nx'*Vx'*Mx*Vx*Nx;

    Ny=null([Btu' Bsu' Dzu']);

    Vy=[     Att'             Ast'            Ctz'; ...
           -eye(nt)       zeros(nt,ns)     zeros(nt,nz); ...
            Ats'             Ass'             Csz'; ...
           zeros(ns,nt)       -eye(ns)       zeros(ns,nz); ...
              Btd'             Bsd'             Dzd'; ...
            zeros(nz,nt)      zeros(nz,ns)       -eye(nz)];
               
    My=blkdiag([zeros(nt) YT{i}; YT{i} zeros(nt)],M2,BlockGammaInv);
    Ly=Ny'*Vy'*My*Vy*Ny;

    Lxy=[XT{i} eye(nt); eye(nt) YT{i}];
    
    LMI=[LMI, Lx<=0, Ly>=0, XT{i}>=0, YT{i}>=0, Lxy>=0, alpha{i}>=0, beta{i}>=0];
    
    Sum=Sum+Sub.Nk*trace(Btd'*XT{i}*Btd);
end

L1=[D eye(ns); eye(ns) Di];
LMI=[LMI,D<=0, Di<=0, L1<=0];
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
    Sub=GP.Sub{i};
    nt=size(Sub.A,1)-ns;
    
    YT{i}=value(YT{i});
    XT{i}=value(XT{i});
    alpha{i}=value(alpha{i});
    beta{i}=value(beta{i});
    
    [Qt,Rt] = qr(eye(nt)-XT{i}*YT{i});
    Xgk=Qt;
    Ygk=Rt';
    Xk=-Xgk'*YT{i}/Ygk';
    Xk=(Xk+Xk')./2;
    XC{i}=[XT{i} Xgk; Xgk' Xk];
end

% Construct closed loop multipliers
M1=value(M1);
M2=value(M2);
M1(isnan(M1))=0;
M2(isnan(M2))=0;

TC=[eye(ns) zeros(ns,3*ns); zeros(ns,2*ns) eye(ns) zeros(ns); ...
    zeros(ns) eye(ns) zeros(ns,2*ns); zeros(ns,3*ns) eye(ns)];

Ma=M1-inv(M2);
M12=Ma;
M22=Ma;
MC=[M1 M12; M12' M22];

MC=TC*MC*TC;
MC11=MC(1:2*ns,1:2*ns);
MC12=MC(1:2*ns,2*ns+1:end);
MC22=MC(2*ns+1:end,2*ns+1:end);

%% Controller Reconstruction

for i=1:ALPHA
    Sub=GP.Sub{i};
    nt=size(Sub.A,1)-ns;
    nu=GP.Sub{i}.nu;
    nd=size(Sub.B,2)-nu;
    ny=GP.Sub{i}.ny;
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
        
    Btuc=[zeros(nt,nt+ns) Btu; eye(nt) zeros(nt,ns+nu)];
    Bsuc=[zeros(ns,nt+ns) Bsu; zeros(ns,nt) eye(ns) zeros(ns,nu)];
    Dzuc=[zeros(nz,nt+ns) Dzu];
    V=[Btuc;Bsuc;Dzuc]';

    Ctyc=[zeros(nt) eye(nt); zeros(ns,2*nt); Cty zeros(ny,nt)];
    Csyc=[zeros(nt,2*ns); zeros(ns)  eye(ns); Csy zeros(ny,ns)];
    Dydc=[zeros(nt+ns,nd); Dyd];
    U=[Ctyc Csyc Dydc];

    R11=[Att zeros(nt) Ats zeros(size(Ats)); zeros(nt,2*nt+2*ns); Ast zeros(ns,nt) Ass zeros(ns,ns); zeros(ns,2*nt+2*ns)];

    R12=[Btd; zeros(size(Btd)); Bsd; zeros(size(Bsd))];
    R21=[Ctz zeros(size(Ctz)) Csz zeros(size(Csz))];
    R22=Dzd;
    R=[R11 R12; R21 R22];
    
    if alpha{i}*beta{i}>=1
        rho{i}=alpha{i};
    else
        rho{i}=1/beta{i};
    end
    
    nc=size(XC{i},1);
    M11=[zeros(nc,nc+2*ns+nd); zeros(2*ns,nc) MC11 zeros(2*ns,nd); zeros(nd,nc+2*ns) -rho{i}*eye(nd)];
    M12=[XC{i} zeros(nc,2*ns+nz); zeros(2*ns,nc) MC12 zeros(2*ns,nz); zeros(nd,nc+2*ns+nz)];
    M22=[zeros(nc,nc+2*ns+nz); zeros(2*ns,nc) MC22 zeros(2*ns,nz); zeros(nz,nc+2*ns) eye(nz)];
    M=[M11 M12; M12' M22];
    
    [~, ~, H]=svd(V);
    [~, ~, J]=svd(U);

    N=blkdiag(J',inv(H))*M*blkdiag(J,inv(H'));
    N=(N+N')./2;
    T=H'*R*J;
    
    n1=nt+ns+nu;
    n2=nt+ns+ny;
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

    n1=nt+ns+nu;
    n2=nt+ns+ny;
    [Vt,Dt] = eig(Ma);
    [~,ind] = sort(diag(Dt));
    Vs =Vt(:,ind);

    Z=Vs(:,1:n2);
    Z1=Z(1:n2,:);
    Z2=Z(n2+1:end,:);
    
    Z11=Z2/Z1;

    V1=V*H;
    V1=V1(:,1:n1);
    U1=U*J;
    U1=U1(:,1:n2);

    Theta=V1'\(Z11-T11)/U1;
    K.Sub{i}.A=Theta(1:nt+ns,1:nt+ns);
    K.Sub{i}.B=Theta(1:nt+ns,nt+ns+1:end);
    K.Sub{i}.C=Theta(nt+ns+1:end,1:nt+ns);
    K.Sub{i}.D=Theta(nt+ns+1:end,nt+ns+1:end);
end

K.Ts=GP.Ts;
K.T=GP.T;
K.Ord=GP.Ord;
K.ns=GP.ns;
gam=sqrt(gammasq);





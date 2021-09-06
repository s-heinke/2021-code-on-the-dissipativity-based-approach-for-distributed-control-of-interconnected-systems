function [K, gam] = h2AG(GP,nmeas,ncon,sub,Meth)
% Distributed H2 control for systems interconnected over an arbitrary graph
% Licensed under GPLv3. See License.txt in the project root for license information.
% Author: Simon Heinke

% Inputs:
% GP: Generalized Plant, with subsystems GP.Sub{i}, graph adjacency matrix
%       GP.Int, a row vector GP.Ord specifying the type of subsystem i and
%       sampling time GP.Ts.
% nmeas: row vector containing the number of measured outputs of subsystem i
% ncon: row vector containing the number of control inputs of subsystem i
% sub: scalar >=1 specifying suboptimality of the bound (sometimes necessary 
%       due to numerical issues)

% Outputs:
% K: interconnected controller
% gam: H2 performance bound

% adjacency matrix
N=GP.Int;
% dimensions
nSub=size(GP.Ord,2);
gammasq=sdpvar(1);

if ~exist('Meth','var')
    Meth=1;
end

% Define the optimization variables
for i=1:nSub
    type=GP.Ord(i);
    ns=sum(GP.Int(:,i));
    Sub=GP.Sub{type};
    nt=size(Sub.A,1)-ns;
    
    YT{i}=sdpvar(nt);
    XT{i}=sdpvar(nt);

    alpha{i}=sdpvar(1);
    beta{i}=sdpvar(1);
    
    for j=1:nSub
        if N(i,j)~=0
            X11{i,j}=sdpvar(N(i,j));
            Y11{i,j}=sdpvar(N(i,j));
            if i>j
               X12{i,j}=sdpvar(N(i,j),N(i,j),'full');
               Y12{i,j}=sdpvar(N(i,j),N(i,j),'full');
            end
        end
    end
end

% Construct the multipliers
for i=1:nSub
    P11{i}=[];
    P12{i}=[];
    P22{i}=[];
    Z11{i}=[];
    Z12{i}=[];
    Z22{i}=[];
    for j=1:nSub
        if N(i,j)~=0
            P11{i}=blkdiag(P11{i},X11{i,j});
            Z11{i}=blkdiag(Z11{i},Y11{i,j});
            if j<i
                P12{i}=blkdiag(P12{i},X12{i,j});
                Z12{i}=blkdiag(Z12{i},Y12{i,j});
            elseif j>i
                P12{i}=blkdiag(P12{i},-X12{j,i}');
                Z12{i}=blkdiag(Z12{i},-Y12{j,i}');
            else
                error('i=j')
            end
            P22{i}=blkdiag(P22{i},-X11{j,i});
            Z22{i}=blkdiag(Z22{i},-Y11{j,i});
        end
    end
end

% Construct the LMIs
LMI=[];
Sum=0;
for i=1:nSub
    type=GP.Ord(i);
    ns=sum(GP.Int(:,i));
    Sub=GP.Sub{type};
    nt=size(Sub.A,1)-ns;
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

    Vx=[      eye(nt)         zeros(nt,ns)     zeros(nt,nd); ...
                Att             Ats              Btd; ...
            zeros(ns,nt)       eye(ns)         zeros(ns,nd); ...
                Ast             Ass             Bsd; ...
            zeros(nd,nt)      zeros(nd,ns)       eye(nd); ...
               Ctz              Csz             Dzd];

    Mx=blkdiag([zeros(nt) XT{i}; XT{i} zeros(nt)],[P11{i} P12{i}; P12{i}' P22{i}],blkdiag(-alpha{i}*eye(nd),eye(nz)));
    Lx=Nx'*Vx'*Mx*Vx*Nx;

    Ny=null([Btu' Bsu' Dzu']);

    Vy=[      Att'             Ast'            Ctz'; ...
            -eye(nt)       zeros(nt,ns)     zeros(nt,nz); ...
              Ats'             Ass'             Csz'; ...
           zeros(ns,nt)       -eye(ns)       zeros(ns,nz); ...
             Btd'             Bsd'             Dzd';...  
            zeros(nz,nt)      zeros(nz,ns)       -eye(nz)];

    My=blkdiag([zeros(nt) YT{i}; YT{i} zeros(nt)],[Z11{i} Z12{i}; Z12{i}' Z22{i}],blkdiag(-beta{i}*eye(nd),eye(nz)));
    Ly=Ny'*Vy'*My*Vy*Ny;

    Lxy=[XT{i} eye(nt); eye(nt) YT{i}];
        
    LMI=[LMI, Lx<=0, Ly>=0, XT{i}>=0, YT{i}>=0, Lxy>=0, alpha{i}>=0, beta{i}>=0];
    
    Sum=Sum+trace(Btd'*XT{i}*Btd);
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

% Construct the closed loop Lyapunov matrix
for i=1:nSub
    type=GP.Ord(i);
    ns=sum(GP.Int(:,i));
    Sub=GP.Sub{type};
    nt=size(Sub.A,1)-ns;
    
    YT{i}=value(YT{i});
    XT{i}=value(XT{i});
    
    [Qt,Rt] = qr(eye(nt)-XT{i}*YT{i});
    XC{i}=[eye(nt) XT{i}; zeros(nt) Qt']/[YT{i} eye(nt); Rt zeros(nt)];
    
    for j=1:nSub
        if N(i,j)~=0
            X11{i,j}=value(X11{i,j});
            Y11{i,j}=value(Y11{i,j});
            if i>j
               X12{i,j}=value(X12{i,j});
               X12{i,j}(isnan(X12{i,j}))=0;
               Y12{i,j}=value(Y12{i,j});
               Y12{i,j}(isnan(Y12{i,j}))=0;
            end
        end
    end
end

% Construct the closed loop multipliers
Case=1; % Case 1: n_{ij}^K=n_{ij}, Case 2: n_{ij}^K=3n_{ij}.
for i=1:nSub
    for j=1:nSub
       if i > j
           if N(i,j)~=0
                S11{i,j}=[X11{i,j} X12{i,j}; X12{i,j}' -X11{j,i}];
                R11{i,j}=[Y11{i,j} Y12{i,j}; Y12{i,j}' -Y11{j,i}];
                T=[S11{i,j} eye(size(S11{i,j},1)); eye(size(S11{i,j},2)) R11{i,j}];
                NoNeg=sum(eig(T)<0); % Number of negative Eigenvalues
                NoPos=size(T,1)-NoNeg;
                if NoNeg~=NoPos
                    Case=2; % if for any egde NoNeg~=NoPos switch to Case 2
                end
           end
       end
    end
end

if Case==1 
    for i=1:nSub
        for j=1:nSub
            if i>j
                if N(i,j)~=0
                    if Meth==1
                        [Q,R]=qr(eye(size(S11{i,j},1))-S11{i,j}*R11{i,j});
                        S12=Q;
                        R12=R';
                        S22=-S12'*R11{i,j}/R12';
                    elseif Meth==2
                        Ma=S11{i,j}-inv(R11{i,j});
                        [V,D]=eig(Ma);
                        [~,I] = sort(diag(D),'descend');
                        V=V(:,I);
                        D=D(I,I);
                        S12=V*sqrt(abs(D));
                        S22=sign(D);
                    elseif Meth==3
                        S12=S11{i,j}-inv(R11{i,j});
                        S22=S12;
                    else
                        error('undefined method')
                    end
                    ns=N(i,j);
                    nsk=ns;
                    
                    X11gk{i,j}=S12(1:ns,1:nsk);
                    X12gk{i,j}=S12(1:ns,nsk+1:end);
                    X12kg{i,j}=S12(ns+1:end,1:nsk)';
                    X11gk{j,i}=-S12(ns+1:end,nsk+1:end);
                    X11k{i,j}=S22(1:nsk,1:nsk);
                    X12k{i,j}=S22(1:nsk,nsk+1:end);
                    X11k{j,i}=-S22(nsk+1:end,nsk+1:end);
                end
            end
        end
    end
elseif Case==2
    for i=1:nSub
        for j=1:nSub
            if i > j
                if N(i,j)~=0
                    T=S11{i,j}-inv(R11{i,j});
                    [V,D]=eig(T);
                    [~,I] = sort(diag(D),'descend');
                    V=V(:,I);
                    D=D(I,I);
                    V=V*sqrt(abs(D));
                    D=sign(D);
                    p=sum(sum(D>0));
                    Vp=V(:,1:p);
                    Vm=V(:,p+1:end);
                    R1p=sum(sign(eig(R11{i,j}))>0);
                    R1m=sum(sign(eig(R11{i,j}))<0);
                    p=4*N(i,j)-R1p;
                    m=4*N(i,j)-R1m;
                    J=blkdiag(eye(p),-eye(m));
                    S13=J;
                    S12=[Vp zeros(size(Vp,1),p-size(Vp,2)) Vm zeros(size(Vm,1),m-size(Vm,2))];
                    ns=N(i,j);
                    nsk=3*ns;
                    X11gk{i,j}=S12(1:ns,1:nsk);
                    X12gk{i,j}=S12(1:ns,nsk+1:end);
                    X12kg{i,j}=S12(ns+1:end,1:nsk)';
                    X11gk{j,i}=-S12(ns+1:end,nsk+1:end);
                    X11k{i,j}=S13(1:nsk,1:nsk);
                    X12k{i,j}=S13(1:nsk,nsk+1:end);
                    X11k{j,i}=-S13(nsk+1:end,nsk+1:end);
                end
            end
        end
    end 
end

for i=1:nSub
    P11G{i}=[];
    P12G{i}=[];
    P22G{i}=[];
    
    P11K{i}=[];
    P12K{i}=[];
    P22K{i}=[];
    
    P11GK{i}=[];
    P12GK{i}=[];
    P22GK{i}=[];
    P12KG{i}=[];
    
    for j=1:nSub
        if N(i,j)~=0
            P11G{i}=blkdiag(P11G{i},X11{i,j});
            P11K{i}=blkdiag(P11K{i},X11k{i,j});
            P11GK{i}=blkdiag(P11GK{i},X11gk{i,j});
            if j<i
                P12G{i}=blkdiag(P12G{i},X12{i,j});
                P12K{i}=blkdiag(P12K{i},X12k{i,j});
                P12GK{i}=blkdiag(P12GK{i},X12gk{i,j});
                P12KG{i}=blkdiag(P12KG{i},X12kg{i,j});
            elseif j>i
                P12G{i}=blkdiag(P12G{i},-X12{j,i}');
                P12K{i}=blkdiag(P12K{i},-X12k{j,i}');
                P12GK{i}=blkdiag(P12GK{i},-X12kg{j,i}');
                P12KG{i}=blkdiag(P12KG{i},-X12gk{j,i}');
            else
                error('i=j')
            end
            P22G{i}=blkdiag(P22G{i},-X11{j,i});
            P22K{i}=blkdiag(P22K{i},-X11k{j,i});
            P22GK{i}=blkdiag(P22GK{i},-X11gk{j,i});
        end
    end
    P11C{i}=[P11G{i} P11GK{i}; P11GK{i}' P11K{i}];
    P12C{i}=[P12G{i} P12GK{i}; P12KG{i} P12K{i}];
    P22C{i}=[P22G{i} P22GK{i}; P22GK{i}' P22K{i}];
end

% Controller Reconstruction
for i=1:nSub
    type=GP.Ord(i);
    ns=sum(GP.Int(:,i));
    Sub=GP.Sub{type};
    nt=size(Sub.A,1)-ns;
    nu=ncon(i);
    nd=size(Sub.B,2)-nu;
    ny=nmeas(i);
    nz=size(Sub.C,1)-ny;
    if Case==1
        nsk=ns;
    elseif Case==2
        nsk=3*ns;
    end
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
        
    Btuc=[zeros(nt,nt+nsk) Btu; eye(nt) zeros(nt,nsk+nu)];
    Bsuc=[zeros(ns,nt+nsk) Bsu; zeros(nsk,nt) eye(nsk) zeros(nsk,nu)];
    Dzuc=[zeros(nz,nt+nsk) Dzu];
    V=[Btuc;Bsuc;Dzuc]';

    Ctyc=[zeros(nt) eye(nt); zeros(nsk,2*nt); Cty zeros(ny,nt)];
    Csyc=[zeros(nt,ns+nsk); zeros(nsk,ns)  eye(nsk); Csy zeros(ny,nsk)];
    Dydc=[zeros(nt+nsk,nd); Dyd];
    U=[Ctyc Csyc Dydc];
    
    R11=[Att zeros(nt) Ats zeros(nt,nsk); zeros(nt,2*nt+ns+nsk); Ast zeros(ns,nt) Ass zeros(ns,nsk); zeros(nsk,2*nt+ns+nsk)];
    R12=[Btd; zeros(size(Btd)); Bsd; zeros(nsk,nd)];
    R21=[Ctz zeros(size(Ctz)) Csz zeros(nz,nsk)];
    R22=Dzd;
    R=[R11 R12; R21 R22];
    
    if value(alpha{i})*value(beta{i})>=1
        rho{i}=value(alpha{i});
    else
        rho{i}=1/value(beta{i});
    end
    nc=size(XC{i},1);
    M11=[zeros(nc,nc+ns+nsk+nd); zeros(ns+nsk,nc) P11C{i} zeros(ns+nsk,nd); zeros(nd,nc+ns+nsk) -rho{i}*eye(nd)];
    M12=[XC{i} zeros(nc,ns+nsk+nz); zeros(ns+nsk,nc) P12C{i} zeros(ns+nsk,nz); zeros(nd,nc+ns+nsk+nz)];
    M22=[zeros(nc,nc+ns+nsk+nz); zeros(ns+nsk,nc) P22C{i} zeros(ns+nsk,nz); zeros(nz,nc+ns+nsk) eye(nz)];
    M=[M11 M12; M12' M22];
    
    [~, ~, H]=svd(V);
    [~, ~, J]=svd(U);

    Ni=blkdiag(J',inv(H))*M*blkdiag(J,inv(H'));
    Ni=(Ni+Ni')./2;
    T=H'*R*J;
    
    n1=nt+nsk+nu;
    n2=nt+nsk+ny;
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

    Ma=Q'*Ni*Q-Q'*Ni*S/(S'*Ni*S)*S'*Ni*Q;
    Ma=(Ma+Ma')./2;

    n1=nt+nsk+nu;
    n2=nt+nsk+ny;
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
    K.Sub{i}.A=Theta(1:nt+nsk,1:nt+nsk);
    K.Sub{i}.B=Theta(1:nt+nsk,nt+nsk+1:end);
    K.Sub{i}.C=Theta(nt+nsk+1:end,1:nt+nsk);
    K.Sub{i}.D=Theta(nt+nsk+1:end,nt+nsk+1:end);
end
if Case==1
    K.Int=GP.Int;
elseif Case==2
    K.Int=3*GP.Int;
end
K.Ord=1:nSub;
K.Ts=0;
gam=sqrt(gammasq);


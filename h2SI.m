function [K, gammasq] = h2SI(GP,nmeas,ncon,sub)
% Distributed H2 control for spatially invariant systems
% Licensed under GPLv3. See License.txt in the project root for license information.
% Author: Simon Heinke

% Inputs:
% GP: Generalized Plant, with GP.blk=[nt np nm], where nt is the temporal 
%       order of the subsystem and np and nm are the spatial order in the 
%       positive and negative direction.
% nmeas: number of measured outputs per subsystem
% ncon: number of control inputs per subsystem
% sub: scalar >=1 specifying suboptimality of the bound (sometimes necessary 
%       due to numerical issues)

% Outputs:
% K: spatially invaraint controller
% gam: H2 performance bound

% Define the optimization variables
Sub=GP;
ns=2*GP.blk(2);
nt=size(Sub.A,1)-ns;
nu=ncon;
nd=size(Sub.B,2)-nu;
ny=nmeas;
nz=size(Sub.C,1)-ny;

YT=sdpvar(nt,nt);
XT=sdpvar(nt,nt);

alpha=sdpvar(1);
beta=sdpvar(1);

nv=ns/2;
X11=sdpvar(nv);
X22=sdpvar(nv);
X12=sdpvar(nv,nv,'full');
Y11=sdpvar(nv);
Y22=sdpvar(nv);
Y12=sdpvar(nv,nv,'full');
gammasq=sdpvar(1);

% Construct the LMIs
LMI=[];    
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

Vx=[   eye(nt)         zeros(nt,ns)     zeros(nt,nd); ...
         Att             Ats              Btd; ...
       zeros(ns,nt)     eye(ns)         zeros(ns,nd); ...  
         Ast              Ass             Bsd; ...
       zeros(nd,nt)      zeros(nd,ns)    eye(nd); ...
          Ctz              Csz             Dzd];

BlockGamma=blkdiag(-alpha*eye(nd),eye(nz));
BlockGammaInv=blkdiag(-beta*eye(nd),eye(nz));

P11=blkdiag(X11,X22);
P12=[zeros(nv) X12; -X12' zeros(nv)];
Mi=[P11 P12; P12' -P11];
Mx=blkdiag([zeros(nt) XT; XT zeros(nt)],Mi,BlockGamma);
Lx=Nx'*Vx'*Mx*Vx*Nx;

Ny=null([Btu' Bsu' Dzu']);

Vy=[      Att'             Ast'            Ctz'; ...
        -eye(nt)       zeros(nt,ns)     zeros(nt,nz); ...
          Ats'             Ass'             Csz'; ...
       zeros(ns,nt)       -eye(ns)       zeros(ns,nz); ...
          Btd'             Bsd'             Dzd' 
       zeros(nz,nt)      zeros(nz,ns)       -eye(nz)];

Z11=blkdiag(Y11,Y22);
Z12=[zeros(nv) Y12; -Y12' zeros(nv)];
Min=[Z11 Z12; Z12' -Z11];
My=blkdiag([zeros(nt) YT; YT zeros(nt)],Min,BlockGammaInv);
Ly=Ny'*Vy'*My*Vy*Ny;

Lxy=[XT eye(nt); eye(nt) YT];

LMI=[LMI, Lx<=0, Ly>=0, XT>=0, YT>=0, Lxy>=0, alpha>=0, beta>=0];

Sum=trace(Btd'*XT*Btd);
 
LMI1=[LMI, Sum<=gammasq];

OPTIONS = sdpsettings('solver', 'sedumi');
OPTIONS.showprogess=1;

% Solve the LMIs:
diagnostics = optimize(LMI1,gammasq,OPTIONS)
gammasq=value(gammasq);
if sub~=1
    % solve the LMIs again with a suboptimal bound due to numerical issues
    gammasq=gammasq*sub;
    LMI2=[LMI, Sum<=gammasq];
    diagnostics = optimize(LMI2,[],OPTIONS)
end

alpha=value(alpha);
beta=value(beta);

% Construct the closed loop Lyapunov matrix
YT=value(YT);
XT=value(XT);

[Qtemp,Rtemp] = qr(eye(nt)-XT*YT);
Xgk=Qtemp;
Ygk=Rtemp';
Xk=-Xgk'*YT/Ygk';
Xk=(Xk+Xk')./2;
XC=[XT Xgk; Xgk' Xk];

% Construct the cloosed loop multipliers
X11=value(X11);
X22=value(X22);
X12=value(X12);
X12(isnan(X12))=0;
S1=[X11 X12; X12' -X22];
Y11=value(Y11);
Y22=value(Y22);
Y12=value(Y12);
Y12(isnan(X12))=0;
R1=[Y11 Y12; Y12' -Y22];
[Qt, Rt]=qr(eye(size(S1,1))-S1*R1);
S2=Qt;
R2=Rt';
S3=-S2'*R1/R2';
T=blkdiag(eye(nv),[zeros(nv) eye(nv); eye(nv) zeros(nv)],eye(nv));
S=T*[S1 S2; S2' S3]*T;
XC11=S(1:ns,1:ns);
XC12=S(1:ns,ns+1:end);
XC22=-S(ns+1:end,ns+1:end);

Z11C=T*blkdiag(XC11,XC22)*T;
Z12C=T*[zeros(ns) XC12; -XC12' zeros(ns)]*T;
Z22C=-Z11C;


%% Controller Reconstruction
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

if alpha*beta>=1
    rho=alpha;
else
    rho=1/beta;
end
nc=size(XC,1);
M11=[zeros(nc,nc+2*ns+nd); zeros(2*ns,nc) Z11C zeros(2*ns,nd); zeros(nd,nc+2*ns) -rho*eye(nd)];
M12=[XC zeros(nc,2*ns+nz); zeros(2*ns,nc) Z12C zeros(2*ns,nz); zeros(nd,nc+2*ns+nz)];
M22=[zeros(nc,nc+2*ns+nz); zeros(2*ns,nc) Z22C zeros(2*ns,nz); zeros(nz,nc+2*ns) eye(nz)];
M=[M11 M12; M12' M22];

[~,~,H]=svd(V);
[~,~,J]=svd(U);

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
K.A=Theta(1:nt+ns,1:nt+ns);
K.B=Theta(1:nt+ns,nt+ns+1:end);
K.C=Theta(nt+ns+1:end,1:nt+ns);
K.D=Theta(nt+ns+1:end,nt+ns+1:end);
K.Ts=0;

K.blk=GP.blk;




function M = SI2MIMO(P_SI, type, Nsub, nmeas, ncon, ML, MR)
% Constructs a MIMO model for a spatially invariant periodic system or finite 
% extend system.

% Inputs:
% P_SI: basic building block of spatially invariant system with
%       P_SI.blk=[nt np nm], where nt is the temporal order of the subsystem
%       and np and nm are the spatial order in the positive and negative
%       direction.
% type: 'finite' or 'cyclic'
% nmeas: dimension of measured output
% ncon: dimension of control input
% ML, MR: left and right boundary matrix

reorder=1;
if Nsub < 2
   error('less than 2 subsystems') 
end

nt=abs(P_SI.blk(1));
np=P_SI.blk(2);
nm=P_SI.blk(3);

if ncon==size(P_SI.B,2)
   reorder=0;
end

nd=size(P_SI.B,2)-ncon;
nz=size(P_SI.C,1)-nmeas;

Att=P_SI.A(1:nt,1:nt);
Ats=P_SI.A(1:nt,nt+1:end);
Bt=P_SI.B(1:nt,:);
Ast=P_SI.A(nt+1:end,1:nt);
Ass=P_SI.A(nt+1:end,nt+1:end);
Bs=P_SI.B(nt+1:end,:);
Ct=P_SI.C(:,1:nt);
Cs=P_SI.C(:,nt+1:end);
Dsub=P_SI.D;
    
ATT=kron(eye(Nsub),Att);
ATS=kron(eye(Nsub),Ats);
AST=kron(eye(Nsub),Ast);
ASS=kron(eye(Nsub),Ass);
BT=kron(eye(Nsub),Bt);
BS=kron(eye(Nsub),Bs);
CT=kron(eye(Nsub),Ct);
CS=kron(eye(Nsub),Cs);
D=kron(eye(Nsub),Dsub);

M=ss(ATT,[ATS BT],[AST; CT], [ASS BS; CS D]);

if nm~=np
    error('error')
end

if type=='finite'
    T1=diag(ones(Nsub-1,1),-1);
    T2=diag(ones(Nsub-1,1),1);
    M1=[eye(nm) zeros(nm); zeros(nm) zeros(nm)];
    M2=[zeros(nm) zeros(nm); zeros(nm) eye(nm)];
    T=kron(T1,M1)+kron(T2,M2);
    T(1:2*nm,1:2*nm)=[zeros(nm) ML; zeros(nm) zeros(nm)];
    T((end-2*nm+1):end,(end-2*nm+1):end)=[zeros(nm) zeros(nm); MR zeros(nm)];
elseif type =='cyclic'
    TL=diag(ones(Nsub-1,1),-1);
    TH=diag(ones(Nsub-1,1),1);
    TH(Nsub,1)=1; TL(1,Nsub)=1;
    ML=[eye(nm) zeros(nm); zeros(nm) zeros(nm)];
    MH=[zeros(nm) zeros(nm); zeros(nm) eye(nm)];
    T=kron(TL,ML)+kron(TH,MH);
end

M=lft(T,M);
% reordering the inputs and outputs in case of generalized plant
if reorder
   T1=kron(eye(Nsub),[eye(nd) zeros(nd,ncon)]);
   T2=kron(eye(Nsub),[zeros(ncon,nd) eye(ncon)]);
   T=[T1; T2];
   M.B=M.B/T;
   M.D=M.D/T;
   T1=kron(eye(Nsub),[eye(nz) zeros(nz,nmeas)]);
   T2=kron(eye(Nsub),[zeros(nmeas,nz) eye(nmeas)]);
   T=[T1; T2];
   M.C=T*M.C;
   M.D=T*M.D;
end

if P_SI.blk(1)>0
    M.Ts=0;
else
    M.Ts=-1;
end


















function G = AS2MIMO(P,ny,nu)
% Construct a MIMO model from an alpha-heterogeneous string interconnected model

% Inputs:
% P: plant with subsystems P.Sub{i}, row vector P.Ord specifying the type
%       of subsystems and P.ns, the dimension of the interconnection signals

reorder=1;
if ~exist('nu','var')
   nu=0;
   reorder=0;
end
if ~exist('ny','var')
   ny=0;
end
nsub=size(P.Ord,2);
ATT=[];
ATS=[];
AST=[];
ASS=[];
BT=[];
BS=[];
CT=[];
CS=[];
D=[];
ns=P.ns;
for i=1:nsub
    if i==1 || i==nsub
        nsi=ns;
    else
        nsi=2*ns;
    end
    type=P.Ord(i);
    Sub=P.Sub{type};
    nt=size(Sub.A,1)-nsi;
    
    Att=Sub.A(1:nt,1:nt);
    Ats=Sub.A(1:nt,nt+1:end);
    Bt=Sub.B(1:nt,:);
    Ast=Sub.A(nt+1:end,1:nt);
    Ass=Sub.A(nt+1:end,nt+1:end);
    Bs=Sub.B(nt+1:end,:);
    Ct=Sub.C(:,1:nt);
    Cs=Sub.C(:,nt+1:end);
    Dsub=Sub.D;
    
    ATT=blkdiag(ATT,Att);
    ATS=blkdiag(ATS,Ats);
    AST=blkdiag(AST,Ast);
    ASS=blkdiag(ASS,Ass);
    BT=blkdiag(BT,Bt);
    BS=blkdiag(BS,Bs);
    CT=blkdiag(CT,Ct);
    CS=blkdiag(CS,Cs);
    D=blkdiag(D,Dsub);
end

M=ss(ATT, [BT ATS], [CT; AST], [D CS; BS ASS],P.Ts);

T=kron(diag(ones(nsub-3,1),-1),[eye(ns) zeros(ns); zeros(ns) zeros(ns)])...
    +kron(diag(ones(nsub-3,1),1),[zeros(ns) zeros(ns); zeros(ns) eye(ns)]);
T=blkdiag(zeros(ns),T,zeros(ns));
T(1:ns,2*ns+1:3*ns)=eye(ns);
T(ns+1:2*ns,1:ns)=eye(ns);
T(end-ns+1:end,end-3*ns+1:end-2*ns)=eye(ns);
T(end-2*ns+1:end-ns,end-ns+1:end)=eye(ns);

G=lft(M,T);
if reorder
   for i=1:nsub
       index=P.Ord(i);
       nd(i)=size(P.Sub{index}.B,2)-nu(i);
       nz(i)=size(P.Sub{index}.C,1)-ny(i);
   end

   in_in=[1:nd(1)]; 
   in_out=[1:nz(1)]; 
   s_nd=0;
   s_nu=0;
   s_nz=0;
   s_ny=0;
   % for z and d 
   for i=2:nsub
       s_nd=s_nd+nd(i-1);
       s_nu=s_nu+nu(i-1);
       in_in=[in_in s_nd+s_nu+1:s_nd+s_nu+nd(i)];
       s_nz=s_nz+nz(i-1);
       s_ny=s_ny+ny(i-1);
       in_out=[in_out s_nz+s_ny+1:s_nz+s_ny+nz(i)];
   end
   % for y and u 
   in_in=[in_in nd(1)+1:nu(1)+nd(1)]; 
   in_out=[in_out nz(1)+1:ny(1)+nz(1)];
   s_nd=0;
   s_nu=0;
   s_nz=0;
   s_ny=0;
   for i=2:nsub
       s_nd=s_nd+nd(i-1);
       s_nu=s_nu+nu(i-1);
       in_in=[in_in s_nd+s_nu+nd(i)+1:s_nd+s_nu+nd(i)+nu(i)];
       s_nz=s_nz+nz(i-1);
       s_ny=s_ny+ny(i-1);
       in_out=[in_out s_nz+s_ny+nz(i)+1:s_nz+s_ny+nz(i)+ny(i)];
   end
   
   G=G(in_out,in_in);
end


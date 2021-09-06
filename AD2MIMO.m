function G = AD2MIMO(P,ny,nu)
% Construct a MIMO model from an alpha-heterogeneous decomposable model

% Inputs:
% P: plant with subsystems GP.Sub{i}, interconnection matrix
%       GP.T, size of the interconnection signals GP.ns and a row vector 
%       GP.Ord specifying the type of subsystem i.
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
    type=P.Ord(i);
    Sub=P.Sub{type};
    nt=size(Sub.A,1)-ns;
    
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

T=P.T;
Ti=kron(T,eye(ns));

G=lft(M,Ti);

if reorder
   for i=1:nsub
       index=P.Ord(i);
       nug(i)=nu(index);
       nyg(i)=ny(index);
       ndg(i)=size(P.Sub{index}.B,2)-nu(index);
       nzg(i)=size(P.Sub{index}.C,1)-ny(index);
   end

   in_in=[1:ndg(1)]; 
   in_out=[1:nzg(1)]; 
   s_nd=0;
   s_nu=0;
   s_nz=0;
   s_ny=0;
   % for z and d 
   for i=2:nsub
       s_nd=s_nd+ndg(i-1);
       s_nu=s_nu+nug(i-1);
       in_in=[in_in s_nd+s_nu+1:s_nd+s_nu+ndg(i)];
       s_nz=s_nz+nzg(i-1);
       s_ny=s_ny+nyg(i-1);
       in_out=[in_out s_nz+s_ny+1:s_nz+s_ny+nzg(i)];
   end
   % for y and u 
   in_in=[in_in ndg(1)+1:nug(1)+ndg(1)]; 
   in_out=[in_out nzg(1)+1:nyg(1)+nzg(1)];
   s_nd=0;
   s_nu=0;
   s_nz=0;
   s_ny=0;
   for i=2:nsub
       s_nd=s_nd+ndg(i-1);
       s_nu=s_nu+nug(i-1);
       in_in=[in_in s_nd+s_nu+ndg(i)+1:s_nd+s_nu+ndg(i)+nug(i)];
       s_nz=s_nz+nzg(i-1);
       s_ny=s_ny+nyg(i-1);
       in_out=[in_out s_nz+s_ny+nzg(i)+1:s_nz+s_ny+nzg(i)+nyg(i)];
   end
   
   G=G(in_out,in_in);
end



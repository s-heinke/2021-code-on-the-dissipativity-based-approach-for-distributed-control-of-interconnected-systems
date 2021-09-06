function G = AG2MIMO(P,nmeas,ncon)
% Construct MIMO model from system interconnected over an arbitrary graph
nsub=size(P.Ord,2);
Ts=P.Ts;
reorder=1;
if ~exist('ncon','var')
   ncon=0;
   reorder=0;
end
if ~exist('nmeas','var')
   nmeas=0;
end
nu=ncon;
ny=nmeas;
ATT=[];
ATS=[];
AST=[];
ASS=[];
BT=[];
BS=[];
CT=[];
CS=[];
D=[];
Nint=zeros(1,nsub);
for i=1:nsub
    type=P.Ord(i);
    ns=sum(P.Int(:,i));
    Sub=P.Sub{type};
    nt=size(Sub.A,1)-ns;
    Nint(i)=ns;
    
    
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

M=ss(ATT, [BT ATS], [CT; AST], [D CS; BS ASS],Ts);

Ti=P.Int;
        
ns=zeros(1,nsub);
for i=1:nsub
    ns(i)=sum(Ti(:,i));
end
nT=sum(ns);
T=zeros(nT); %initialize T

for i=1:nsub
    for j=i+1:nsub
        if Ti(j,i)>0
            nj=sum(ns(1:j-1))+1; %location of Tji in T
            ni=sum(ns(1:i-1))+1;
            addj=sum(Ti(j,1:i-1));
            addi=sum(Ti(i,1:j-1));
            T(nj+addj:nj+addj+Ti(j,i)-1,ni+addi:ni+addi+Ti(j,i)-1)=eye(Ti(j,i));
            T(ni+addi:ni+addi+Ti(j,i)-1,nj+addj:nj+addj+Ti(j,i)-1)=eye(Ti(j,i));
        end
    end
end     

nv=size(T,2);
G=lft(M,T,nv,nv);

% reorder inputs and outputs
if reorder
   for i=1:nsub
       index=P.Ord(i);
       nd(i)=size(P.Sub{1,index}.B,2)-nu(i);
       nz(i)=size(P.Sub{1,index}.C,1)-ny(i);
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





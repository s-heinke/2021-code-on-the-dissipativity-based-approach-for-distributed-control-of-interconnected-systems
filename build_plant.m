function [Gp, Gp_SI, Gp_AG, Gp_AD, Gp_AD2, Gp_de] = build_plant(p,N)
    % This file constructs the generalized plants for the interconnected 
    % seesaw-cart system.
    % Author: Simon Heinke
    
    % cost function
    Dzu=0.5;
    Ctz=[1 0 0 0];

    %% 1) Spatially Invariant
    Att=[0 0 1 0; ...
        0 0 0 1; ...
          p.a*p.mA*p.g/p.J-(p.k*(p.L^2)/(2*p.J)) -p.mC*p.g/p.J -p.cappa/p.J -p.c*p.Ke*p.Kt/(p.J*p.Ra*(p.p^2)*(p.r^2)); ...
          p.c*p.a*p.mA*p.g/p.J-p.g-p.c*p.k*p.L^2/(2*p.J) -p.c*p.mC*p.g/p.J -p.c*p.cappa/p.J -(p.c^2/p.J+1/p.mC)*p.Ke*p.Kt/(p.Ra*p.p^2*p.r^2)];
      
    Ats=[0 0; 0 0; -p.k*(p.L^2)/(4*p.J) -p.k*(p.L^2)/(4*p.J); -p.c*p.k*p.L^2/(4*p.J) -p.c*p.k*p.L^2/(4*p.J)];
    Ast=[1 0 0 0; 1 0 0 0];
    Ass=zeros(2);
    Btd=-[0; 0; 1/p.J; p.c/p.J];
    Bsd=zeros(2,1);
    Btu=[0; 0; p.c*p.Kt/(p.J*p.Ra*p.p*p.r); (p.c^2/p.J+1/p.mC)*p.Kt/(p.Ra*p.p*p.r)];
    Bsu=zeros(2,1);
    Csz=zeros(1,2);
    Cty=[eye(2) zeros(2)];        
    Csy=zeros(2);
    Dyd=zeros(2,1);
    Dyu=zeros(2,1);
    Dzd=0;

    % Generalized Plant
    Gp_SI.A=[Att Ats; Ast Ass];
    Gp_SI.B=[Btd Btu; Bsd Bsu];
    Gp_SI.C=[Ctz Csz; Cty Csy];
    Gp_SI.D=[Dzd Dzu; Dyd Dyu];
    Gp_SI.blk=[4 1 1];
    
    %% 2) Arbitrary Graph
    % Subsystem 1 and 3 are identical, subsystem 2 is identical to the
    % subsystem in the spatially invariant case
    Att1=[0 0 1 0; ...
        0 0 0 1; ...
          p.a*p.mA*p.g/p.J-(p.k*(p.L^2)/(4*p.J)) -p.mC*p.g/p.J -p.cappa/p.J -p.c*p.Ke*p.Kt/(p.J*p.Ra*(p.p^2)*(p.r^2)); ...
          p.c*p.a*p.mA*p.g/p.J-p.g-p.c*p.k*p.L^2/(4*p.J) -p.c*p.mC*p.g/p.J -p.c*p.cappa/p.J -(p.c^2/p.J+1/p.mC)*p.Ke*p.Kt/(p.Ra*p.p^2*p.r^2)];
    Ats1=[0; 0; -p.k*(p.L^2)/(4*p.J); -p.c*p.k*p.L^2/(4*p.J)];
    Ast1=[1 0 0 0];
    Ass1=0;
    Btd1=Btd;
    Bsd1=0;
    Btu1=Btu;
    Bsu1=0;
    Ctz1=Ctz;
    Csz1=0;
    Dzd1=Dzd;
    Dzu1=Dzu;
    Cty1=Cty;
    Csy1=zeros(2,1);
    Dyd1=Dyd;
    Dyu1=Dyu;

    % Generalized Plant
    G1.A=[Att1 Ats1; Ast1 Ass1];
    G1.B=[Btd1 Btu1; Bsd1 Bsu1];
    G1.C=[Ctz1 Csz1; Cty1 Csy1];
    G1.D=[Dzd1 Dzu1; Dyd1 Dyu1];                     

    G2.A=[Att Ats; Ast Ass];
    G2.B=[Btd Btu; Bsd Bsu];
    G2.C=[Ctz Csz; Cty Csy];
    G2.D=[Dzd Dzu; Dyd Dyu];                         

    Gp_AG.Sub={G1, G2};                
    Gp_AG.Int=diag(ones(N-1,1),1)+diag(ones(N-1,1),-1);  
    Ord=2*ones(1,N); Ord(1)=1; Ord(N)=1;
    Gp_AG.Ord=Ord;                 
    Gp_AG.Ts=0;
    
    %% 3 a) alpha-heterogeneous decomposable
    G1.nu=1;
    G1.ny=2;
    G1.Nk=3;
    G1.A=[Att1 2*Ats1; Ast1 2*Ass1];
    Gp_AD.Sub={G1};
    P=eye(N)+diag(ones(N-1,1),1)+diag(ones(N-1,1),-1);
    P(1,1)=0; P(N,N)=0;
    Gp_AD.T=0.5*P;
    Gp_AD.ns=1;
    Gp_AD.Ord=ones(1,N);
    Gp_AD.Ts=0;
    %% 3 b) alpha-heterogeneous decomposable
    G1.nu=1;
    G1.ny=2;
    G1.Nk=2;
    G3=G1;
    G3.nu=1;
    G3.ny=2;
    G3.Nk=1;
    P=diag(ones(N-1,1),1)+diag(ones(N-1,1),-1);
    scale=max(abs(eig(P)));
    G1.A=[Att1 scale*Ats1; Ast1 scale*Ass1];
    G3.A=[Att scale*Ats1; Ast1 scale*Ass1];
    Gp_AD2.Sub={G1,G3};
    Gp_AD2.T=(1/scale)*P;
    Gp_AD2.ns=1;
    Gp_AD2.Ord=[1 2 1];
    Gp_AD2.Ts=0;
    %% 4) centralized
    nmeas=2;       
    ncon=1;
    ML=-1;
    MR=-1;
    type='finite';
    % Constructing generalized
    Gp = SI2MIMO(Gp_SI, type, N, nmeas, ncon, ML, MR);
    %% 5) Decentralized
    A=[0 0 1 0; ...
       0 0 0 1; ...
       p.a*p.mA*p.g/p.J    -p.mC*p.g/p.J    -p.cappa/p.J    -p.c*p.Ke*p.Kt/(p.J*p.Ra*(p.p^2)*(p.r^2)); ...
       p.c*p.a*p.mA*p.g/p.J-p.g    -p.c*p.mC*p.g/p.J   -p.c*p.cappa/p.J    -(p.c^2/p.J+1/p.mC)*p.Ke*p.Kt/(p.Ra*p.p^2*p.r^2)];

    Gp_de=ss(A,[Btd Btu],[Ctz; Cty],[Dzd Dzu; Dyd Dyu]);
    
end



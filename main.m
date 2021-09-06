%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% For paper, "On the Dissipativity Based Approach for Distributed Control of Interconnected Systems" by S. Heinke and H. Werner.
% Copyright (c) Institute of Control Systems, Hamburg University of Technology. All rights reserved.
% Licensed under GPLv3. See License.txt in the project root for license information.
% Author: Simon Heinke
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all; close all; clc

% Load parameters
load('parameters')
N=3;    % No. of subsystems
nui=1;
nyi=2;

% Construct generalized plant
% z_i=rho_i+0.5u_i 
[Gp, Gp_SI, Gp_AG, Gp_AD, Gp_AD2, Gp_de]=build_plant(p,N);
%% Controller Synthesis
sub=1.05;
% Centralized Controller
K0=h2syn(Gp,N*nyi,N*nui);
% decentralized controller
K1s=h2syn(Gp_de,nyi,nui);
K1=blkdiag(K1s,K1s,K1s);
% Systems Interconnected over Arbitrary Graphs
ny=ones(1,N)*nyi;
nu=ones(1,N)*nui;
[K_AG, gam1]= h2AG(Gp_AG,ny,nu,sub);
K2 = AG2MIMO(K_AG);

%% Spatially Reversible
[K_SI, gammasq] = h2SI(Gp_SI,nyi,nui,sub);
gam2=sqrt(N*gammasq);
% Construct a spatially reversible controller
K3 = SI2MIMO(K_SI, 'finite', N, nui, nyi,-1, -1);

%% alpha-heterogeneous decomposable
[K_AD, gam3]=h2AD(Gp_AD,sub);
K4=AD2MIMO(K_AD);

[K_AD2, gam4]=h2AD(Gp_AD2,sub);
K5=AD2MIMO(K_AD2);

%% alpha-heterogeneous string interconnected systems
Gp_AS=Gp_AG;
Gp_AS.Sub{3}=Gp_AS.Sub{1};
Gp_AS.Ord(3)=3;
Gp_AS.ns=1;
Gp_AS.Nki=ones(1,3);
[K_AS, gam5]=h2AS(Gp_AS,2*ones(1,3),ones(1,3),sub);
K6=AS2MIMO(K_AS);

%% Different Controller Reconstructions
[K_AG2, gam6]= h2AG(Gp_AG,ny,nu,sub,2);
[K_AG3, gam7]= h2AG(Gp_AG,ny,nu,sub,3);
K7 = AG2MIMO(K_AG2);
K8 = AG2MIMO(K_AG3);
%% Comparison
CL0 = lft(Gp, K0); 
disp(['H2 performance h2syn: ', num2str(norm(CL0,2),3)])

CL0 = lft(Gp, K1); 
disp(['H2 performance h2syn_d: ', num2str(norm(CL0,2),4)])

CL0 = lft(Gp, K2); 
disp(['H2 performance h2AG: ', num2str(norm(CL0,2),3), '         Bound: ', num2str(gam1,4)])

CL0 = lft(Gp, K3); 
disp(['H2 performance h2SR: ', num2str(norm(CL0,2),3), '       Bound: ', num2str(gam2,4)])

CL0 = lft(Gp, K4); 
disp(['H2 performance h2AD1: ', num2str(norm(CL0,2),4), '     Bound: ', num2str(gam3,4)])

CL0 = lft(Gp, K5); 
disp(['H2 performance h2AD2: ', num2str(norm(CL0,2),3), '     Bound: ', num2str(gam4,4)])

CL0 = lft(Gp, K6); 
disp(['H2 performance h2AS: ', num2str(norm(CL0,2),3), '        Bound: ', num2str(gam5,4)])

%% Comparison between different controller reconstructions
CL0 = lft(Gp, K2); 
disp(['H2 performance h2AG (Meth=1): ', num2str(norm(CL0,2),3), '       Bound: ', num2str(gam1,4)])

CL0 = lft(Gp, K7); 
disp(['H2 performance h2AG (Meth=2): ', num2str(norm(CL0,2),3), '     Bound: ', num2str(gam6,4)])

CL0 = lft(Gp, K8); 
disp(['H2 performance h2AG (Meth=3): ', num2str(norm(CL0,2),3), '     Bound: ', num2str(gam7,4)])

%% Comparison in time domain
% disturbance of [0.5 -1 0.2] Nm at time t=1s and duration of 0.5s.
K=K0;
Kd=1;
t_stop=8;
sim('closed_loop',t_stop);
load u_sim
load y_sim
y_limit=[-20 20];
u_limit=[-5 5];
figure(1)
subplot(2,7,1)
    plot(y_sim(1,:),(180/pi)*y_sim(2,:),'b')
    title('h2syn')
    hold on
    grid on
    plot(y_sim(1,:),(180/pi)*y_sim(4,:),'r')
    plot(y_sim(1,:),(180/pi)*y_sim(6,:),'color',[0 0.7 0])
    legend('\rho_1','\rho_2','\rho_3')
    ylabel('angle [°]')
    ylim(y_limit)
    xlim([0 t_stop])
subplot(2,7,8)
    plot(u_sim(1,:),u_sim(2,:),'b')
    title('h2syn')
    hold on
    grid on
    plot(u_sim(1,:),u_sim(3,:),'r')
    plot(u_sim(1,:),u_sim(4,:),'color',[0 0.7 0])
    legend('u_1','u_2','u_3')
    ylabel('voltage [V]')
    xlabel('time [s]')
    ylim(u_limit)
    xlim([0 t_stop])
K=K1;
sim('closed_loop',t_stop);
load u_sim
load y_sim
figure(1)
subplot(2,7,2)
    plot(y_sim(1,:),(180/pi)*y_sim(2,:),'b')
    title('h2syn_d')
    hold on
    grid on
    plot(y_sim(1,:),(180/pi)*y_sim(4,:),'r')
    plot(y_sim(1,:),(180/pi)*y_sim(6,:),'color',[0 0.7 0])
    ylim(y_limit)
    xlim([0 t_stop])
subplot(2,7,9)
    plot(u_sim(1,:),u_sim(2,:),'b')
    title('h2syn_d')
    hold on
    grid on
    plot(u_sim(1,:),u_sim(3,:),'r')
    plot(u_sim(1,:),u_sim(4,:),'color',[0 0.7 0])
    xlabel('time [s]')
    ylim(u_limit)
    xlim([0 t_stop])
K=K2;
sim('closed_loop',t_stop);
load u_sim
load y_sim
figure(1)
subplot(2,7,3)
    plot(y_sim(1,:),(180/pi)*y_sim(2,:),'b')
    title('AG')
    hold on
    grid on
    plot(y_sim(1,:),(180/pi)*y_sim(4,:),'r')
    plot(y_sim(1,:),(180/pi)*y_sim(6,:),'color',[0 0.7 0])
    ylim(y_limit)
    xlim([0 t_stop])
subplot(2,7,10)
    plot(u_sim(1,:),u_sim(2,:),'b')
    title('AG')
    hold on
    grid on
    plot(u_sim(1,:),u_sim(3,:),'r')
    plot(u_sim(1,:),u_sim(4,:),'color',[0 0.7 0])
    xlabel('time [s]')
    ylim(u_limit)
    xlim([0 t_stop])
K=K3;
sim('closed_loop',t_stop);
load u_sim
load y_sim
figure(1)
subplot(2,7,4)
    plot(y_sim(1,:),(180/pi)*y_sim(2,:),'b')
    title('SR')
    hold on
    grid on
    plot(y_sim(1,:),(180/pi)*y_sim(4,:),'r')
    plot(y_sim(1,:),(180/pi)*y_sim(6,:),'color',[0 0.7 0])
    ylim(y_limit)
    xlim([0 t_stop])
subplot(2,7,11)
    plot(u_sim(1,:),u_sim(2,:),'b')
    title('SR')
    hold on
    grid on
    plot(u_sim(1,:),u_sim(3,:),'r')
    plot(u_sim(1,:),u_sim(4,:),'color',[0 0.7 0])
    xlabel('time [s]')
    ylim(u_limit)
    xlim([0 t_stop])
    
K=K4;
sim('closed_loop',t_stop);
load u_sim
load y_sim
figure(1)
subplot(2,7,5)
    plot(y_sim(1,:),(180/pi)*y_sim(2,:),'b')
    title('AD_1')
    hold on
    grid on
    plot(y_sim(1,:),(180/pi)*y_sim(4,:),'r')
    plot(y_sim(1,:),(180/pi)*y_sim(6,:),'color',[0 0.7 0])
    ylim(y_limit)
    xlim([0 t_stop])
subplot(2,7,12)
    plot(u_sim(1,:),u_sim(2,:),'b')
    title('AD_1')
    hold on
    grid on
    plot(u_sim(1,:),u_sim(3,:),'r')
    plot(u_sim(1,:),u_sim(4,:),'color',[0 0.7 0])
    xlabel('time [s]')
    ylim(u_limit)
    xlim([0 t_stop])
    
K=K5;
sim('closed_loop',t_stop);
load u_sim
load y_sim
figure(1)
subplot(2,7,6)
    plot(y_sim(1,:),(180/pi)*y_sim(2,:),'b')
    title('AD_2')
    hold on
    grid on
    plot(y_sim(1,:),(180/pi)*y_sim(4,:),'r')
    plot(y_sim(1,:),(180/pi)*y_sim(6,:),'color',[0 0.7 0])
    ylim(y_limit)
    xlim([0 t_stop])
subplot(2,7,13)
    plot(u_sim(1,:),u_sim(2,:),'b')
    title('AD_2')
    hold on
    grid on
    plot(u_sim(1,:),u_sim(3,:),'r')
    plot(u_sim(1,:),u_sim(4,:),'color',[0 0.7 0])
    xlabel('time [s]')
    ylim(u_limit)
    xlim([0 t_stop])
K=K6;
sim('closed_loop',t_stop);
load u_sim
load y_sim
figure(1)
subplot(2,7,7)
    plot(y_sim(1,:),(180/pi)*y_sim(2,:),'b')
    title('AS')
    hold on
    grid on
    plot(y_sim(1,:),(180/pi)*y_sim(4,:),'r')
    plot(y_sim(1,:),(180/pi)*y_sim(6,:),'color',[0 0.7 0])
    ylim(y_limit)
    xlim([0 t_stop])
subplot(2,7,14)
    plot(u_sim(1,:),u_sim(2,:),'b')
    title('AS')
    hold on
    grid on
    plot(u_sim(1,:),u_sim(3,:),'r')
    plot(u_sim(1,:),u_sim(4,:),'color',[0 0.7 0])
    xlabel('time [s]')
    ylim(u_limit)
    xlim([0 t_stop])
   
K=K7;
sim('closed_loop',t_stop);
load u_sim
load y_sim
figure(2)
subplot(2,2,1)
    plot(y_sim(1,:),(180/pi)*y_sim(2,:),'b')
    title('AG_2')
    hold on
    grid on
    plot(y_sim(1,:),(180/pi)*y_sim(4,:),'r')
    plot(y_sim(1,:),(180/pi)*y_sim(6,:),'color',[0 0.7 0])
    legend('\rho_1','\rho_2','\rho_3')
    ylabel('angle [°]')
    ylim(y_limit)
    xlim([0 t_stop])
subplot(2,2,3)
    plot(u_sim(1,:),u_sim(2,:),'b')
    title('AG_2')
    hold on
    grid on
    plot(u_sim(1,:),u_sim(3,:),'r')
    plot(u_sim(1,:),u_sim(4,:),'color',[0 0.7 0])
    legend('u_1','u_2','u_3')
    ylabel('voltage [V]')
    xlabel('time [s]')
    ylim(u_limit)
    xlim([0 t_stop])
    
K=K8;
sim('closed_loop',t_stop);
load u_sim
load y_sim
figure(2)
subplot(2,2,2)
    plot(y_sim(1,:),(180/pi)*y_sim(2,:),'b')
    title('AG_3')
    hold on
    grid on
    plot(y_sim(1,:),(180/pi)*y_sim(4,:),'r')
    plot(y_sim(1,:),(180/pi)*y_sim(6,:),'color',[0 0.7 0])
    legend('\rho_1','\rho_2','\rho_3')
    ylabel('angle [°]')
    ylim(y_limit)
    xlim([0 t_stop])
subplot(2,2,4)
    plot(u_sim(1,:),u_sim(2,:),'b')
    title('AG_3')
    hold on
    grid on
    plot(u_sim(1,:),u_sim(3,:),'r')
    plot(u_sim(1,:),u_sim(4,:),'color',[0 0.7 0])
    legend('u_1','u_2','u_3')
    ylabel('voltage [V]')
    xlabel('time [s]')
    ylim(u_limit)
    xlim([0 t_stop])
    
    
% set the parameters of the seesaw cart system
clear all; close all; clc

% Physical parameters of the System
a=0.116; c=0.18; L=1; d=0.095; cappa=0.1; p=1/3.71;
Ke=0.0077; L0=0.1; Ls=0.25; k=14;
mA=4.85; mC=1; J=0.68; g=9.81; r=0.0077; Ra=2.6; Kt=0.0077;

para.a=a; para.c=c; para.L=L; para.d=d; para.cappa=cappa; para.p=p;
para.Ke=Ke; para.L0=L0; para.Ls=Ls; para.k=k; para.mA=mA; para.mC=mC;
para.J=J; para.g=g; para.r=r; para.Ra=Ra; para.Kt=Kt;

p=para;
save('parameters.mat','p');
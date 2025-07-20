function [Z] = SCC_Impedancia_Aprox(scc,s,rho)

% UNTITLED Summary of this function goes here
% Detailed explanation goes here

%% ELEMENTOS DA MATRIZ DE IMPEDÂNCIA LONGITUDINAL 
uo=pi*4E-7; jwu=s*uo; 
m.cor=sqrt(jwu/rho.c); 
m.sht=sqrt(jwu/rho.s);
m.arm=sqrt(jwu/rho.a); 

% Sheath thickness
tck.sht=scc.c-scc.b; 

% Armor thickness
tck.arm=scc.e-scc.d;
     
% Impedância interna do núcleo
Z.c = rho.c*m.cor/2/pi/scc.a*coth(0.777*m.cor*scc.a)+0.356*rho.c/pi/(scc.a^2);

% Impedância da blindagem para retorno em sua face interna
Z.si = rho.s*m.sht/2/pi/scc.b*coth(m.sht*tck.sht)-rho.s/2/pi/scc.b/(scc.b+scc.c);

% Impedância da blindagem para retorno em sua face externa
Z.so = rho.s*m.sht/2/pi/scc.c*coth(m.sht*tck.sht)+rho.s/2/pi/scc.c/(scc.b+scc.c);

% IMPEDÂNCIA DE TRANSFERÊNCIA DA BLINDAGEM
% Representa a queda de tensão na face externa da blindagem causada pela 
% circulação de corrente em sua face interna e vice-versa.
% Z12 = Z21 = -Zshth_mut [Pg10-Cp5]
Z.sm = rho.s*m.sht/pi/(scc.b+scc.c)*csch(m.sht*tck.sht);

% Impedância da armadura para retorno em sua face interna
Z.ai = rho.a*m.arm/2/pi/scc.d*coth(m.arm*tck.arm)-rho.a/2/pi/scc.d/(scc.d+scc.e);

% Impedância da armadura para retorno em sua face externa
Z.ao = rho.a*m.arm/2/pi/scc.e*coth(m.arm*tck.arm)+rho.a/2/pi/scc.e/(scc.d+scc.e);

% IMPEDÂNCIA DE TRANSFERÊNCIA DA ARMADURA
% Representa a queda de tensão na face externa da armadura causada pela 
% circulação de corrente em sua face interna e vice-versa.
% Z23 = Z32 = -Zarm_mut [Pg11-Cp5]
Z.am = rho.a*m.arm/pi/(scc.d+scc.e)*csch(m.arm*tck.arm);

% Impedância da isolação primária entre core - sheath
Z.ip = jwu/2/pi*log(scc.b/scc.a);

% Impedância da isolação secundária entre sheath - armor
Z.is = jwu/2/pi*log(scc.d/scc.c);

% Impedância da isolação externa entre armor - soil
Z.ie = jwu/2/pi*log(scc.f/scc.e);
end
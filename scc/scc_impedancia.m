function [Z] = SCC_Impedancia(scc, s, rho)

% UNTITLED Summary of this function goes here
% Detailed explanation goes here

%% CONSTANTES FÍSICAS
uo=pi*4E-7;           % Permeabilidade magnética vácuo
jwu=s*uo;   	        % Frequência angular

%% ELEMENTOS DA MATRIZ DE IMPEDÂNCIA LONGITUDINAL 
m.c=sqrt(jwu/rho.c); m.s=sqrt(jwu/rho.s); m.a=sqrt(jwu/rho.a);        
     
% Impedância interna do núcleo
Z.c = rho.c*m.c/2/pi/scc.a*besseli(0,m.c*scc.a)/besseli(1,m.c*scc.a);

% Impedância da isolação primária entre núcleo e blindagem
Z.ip = jwu/2/pi*log(scc.b/scc.a);

% Impedância da blindagem para retorno em sua face interna
N1 = besseli(0,m.s*scc.b)*besselk(1,m.s*scc.c) + besseli(1,m.s*scc.c)*besselk(0,m.s*scc.b);
D1 = besselk(1,m.s*scc.b)*besseli(1,m.s*scc.c) - besseli(1,m.s*scc.b)*besselk(1,m.s*scc.c);
Z.si = rho.s*m.s/2/pi/scc.b*N1/D1;

% Impedância da blindagem para retorno em sua face externa
N2 = besseli(0,m.s*scc.c)*besselk(1,m.s*scc.b) + besseli(1,m.s*scc.b)*besselk(0,m.s*scc.c);
Z.so = rho.s*m.s/2/pi/scc.c*N2/D1;

% Impedância da isolação secundária entre blindagem e armadura
Z.is = jwu/2/pi*log(scc.d/scc.c);

% Impedância da armadura para retorno em sua face interna
N3 = besseli(0,m.a*scc.d)*besselk(1,m.a*scc.e) + besseli(1,m.a*scc.e)*besselk(0,m.a*scc.d);
D3 = besselk(1,m.a*scc.d)*besseli(1,m.a*scc.e) - besseli(1,m.a*scc.d)*besselk(1,m.a*scc.e);
Z.ai = rho.a*m.a/2/pi/scc.d*N3/D3;

% Impedância da armadura para retorno em sua face externa
N4 = besseli(0,m.a*scc.e)*besselk(1,m.a*scc.d) + besseli(1,m.a*scc.d)*besselk(0,m.a*scc.e);
Z.ao = rho.a*m.a/(2*pi*scc.e)*N4/D3;

% Impedância da isolação externa entre armadura e o solo
Z.ie = jwu/2/pi*log(scc.f/scc.e);
 
%% IMPEDÂNCIA DE TRANSFERÊNCIA DA BLINDAGEM
% Representa a queda de tensão na face externa da blindagem causada pela 
% circulação de corrente em sua face interna e vice-versa.
% Z12 = Z21 = -Zshth_mut [Pg10-Cp5]
Z.sm = rho.s/2/pi/scc.b/scc.c/...
                   	(besselk(1,m.s*scc.b)*besseli(1,m.s*scc.c) - ...
                     besseli(1,m.s*scc.b)*besselk(1,m.s*scc.c));

%% IMPEDÂNCIA DE TRANSFERÊNCIA DA ARMADURA
% Representa a queda de tensão na face externa da armadura causada pela 
% circulação de corrente em sua face interna e vice-versa.
% Z23 = Z32 = -Zarm_mut [Pg11-Cp5]
Z.am = rho.a/2/pi/scc.d/scc.e/...
                    (besselk(1,m.a*scc.d)*besseli(1,m.a*scc.e) - ...
                     besseli(1,m.a*scc.d)*besselk(1,m.a*scc.e));
end
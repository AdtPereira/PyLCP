function [Z,Y] = OHTL_Parametros_Mono(h,ro,sigc,rhog,erg,s)
% Parâmetros de Entrada
% h: altura, m
% ro: raio do condutor, m
% sigc: condutividade do condutor, S/m
% rhog: resistividade do solo, ohms.m
% erg: permissividade relativa do solo,
% s: frequência complexa (c+iw)

% Parâmetros de Saída
% Z: impedância longitudinal, ohms/m
% Y: admitância transversal, S/m

mu0=4*pi*1e-7;
e0=8.854e-12; 
gamag=sqrt(mu0*s.*(1/rhog+s*erg*e0));
p=1./gamag;

Zi=1/(sigc*pi*ro^2)+1/(2*pi*ro)*sqrt(mu0*s/sigc);
Ze=s*mu0/(2*pi)*log(2*h/ro);
Zg=s*mu0/(2*pi).*log((h+p)/h);

Ye=s*2*pi*e0/log(2*h/ro);
Yg=gamag.^2./Zg;

Z=Zi+Ze+Zg;
Y=1./(1./Ye+1./Yg);
end


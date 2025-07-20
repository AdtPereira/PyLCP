%% DADOS DE ENTRADA DA LT

 clc; clear;
   h = 10;                  % Altura [m]
  ri = 0;                   % Raio interno do condutor [m] 
  ro = 0.01;                % Raio externo do condutor [m]
sigc = 1/1.68E-8;           % Condutividade do condutor [S/m]
   f = logspace(3,10,300);  % Faixa de frequência
  jw = 1j*2*pi*f;           % Frequência angular
    
%% Parametros 
FLAG = 3;
APLOG = 0;
[Za,Ya,Zca,Yca,Gammaa] = ohtl_Parametros_QuaseTEM(jw,h,ro,ri,sigc,100,1,FLAG,APLOG);                   
[Zb,Yb,Zcb,Ycb,Gammab] = ohtl_Parametros_QuaseTEM(jw,h,ro,ri,sigc,100,20,FLAG,APLOG); 
[Zc,Yc,Zcc,Ycc,Gammac] = ohtl_Parametros_QuaseTEM(jw,h,ro,ri,sigc,2000,1,FLAG,APLOG); 

%% FIGURA 5.1.A - GRÁFICO DA CONSTANTE DE ATENUAÇÃO EM FUNÇÃO DA FREQUÊNCIA 

figure;
loglog(f,real(squeeze(Za))*1E3,'k',...
       f,real(squeeze(Zb))*1E3,'--k',...
       f,real(squeeze(Zc))*1E3,'-.k');
xlim([1E3 1E10]);
ylim([1E0 1E5]);
grid("on");


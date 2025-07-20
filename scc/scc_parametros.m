function [Zs,Pg,Zg,Yt] = SCC_Parametros(n,s,r,p,rho,er,FlagZg)
% -------------------------------------------------------------------------
% (C) ADILTON JUNIO LADEIRA PEREIRA - PPGEE/UFMG
% (C) ALBERTO RESENDE DE CONTI - PPGEE/UFMG
%
% BELO HORIZONTE, 18 DE DEZEMBRO DE 2023
% -------------------------------------------------------------------------
% Rotina para o cálculo de parâmetros de cabos subterrâneos do tipo single  
% core (SCC). 

%% MATRIZ DE IMPEDÂNCIA LONGITUDINAL
FlagZgAprox=1;
if FlagZgAprox
    [z] = SCC_Impedancia_Aprox(r,s,rho);
else
    [z] = SCC_Impedancia(r,s,rho);
end
 
%% MATRIZ DE IMPEDÂNCIA DE RETORNO PELO SOLO
[Zg] = SCC_Impedancia_Zg(n,p,s,rho.g,er.g,FlagZg);

%% MATRIZ DE IMPEDÂNCIA EM LOOP - ZL (NODA,2008)
if n.cc==1
    zL = z.c+z.ip+Zg(1,1);    
elseif n.cc==2
    zL = [z.c+z.ip+z.si        -z.sm; 
                  -z.sm         z.so+z.ie+Zg(1,1)]; 
elseif n.cc==3
    zL = [z.c+z.ip+z.si        -z.sm                                0; 
                  -z.sm         z.so+z.is+z.ai                  -z.am;
                      0        -z.am               z.ao+z.ie+Zg(1,1)];
end
    
%% ELEMENTOS MÚTUOS DA MATRIZ IMPEDÂNCIA
% MATRIZ DE TRANSFORMAÇÃO A (NODA,2008)
a = eye(n.cc) + diag(-1*ones(n.cc-1,1),-1);

% Forulação principal
if n.f==1 
    A=a; ZL=zL;

elseif n.f==2    
    A =  [          a zeros(n.cc); 
          zeros(n.cc)          a]; 
    
    zL12=zeros(n.cc); zL12(end)=Zg(1,2); 

    ZL = [  zL zL12;
          zL12  zL];
      
elseif n.f==3
    A = [          a        zeros(n.cc)     zeros(n.cc);
         zeros(n.cc)                  a     zeros(n.cc); 
         zeros(n.cc)        zeros(n.cc)              a];
         
   zL12=zeros(n.cc); zL12(end)=Zg(1,2); 
   zL13=zeros(n.cc); zL13(end)=Zg(1,3);
   zL23=zeros(n.cc); zL23(end)=Zg(2,3);
   
   ZL = [  zL	zL12	zL13;
         zL12	  zL	zL23;
         zL13   zL23	 zL];
end

%% ELEMENTOS DA MATRIZ DE IMPEDÂNCIA SÉRIE - Zs (NODA,2008)
Zs = A.'\ZL/A;

%% ADMITÂNCIA TRANSVERSAL 

% Despreza-se o efeito do solo
[ye] = SCC_Admitancia_BaixasFreq(r,s,er,n);

% Matriz Y diagonal em blocos
if n.f==1
    Ye=ye; 

elseif n.f==2    
    Ye = [           ye zeros(n.cc); 
            zeros(n.cc)         ye];
  
elseif n.f==3
    
    Ye = [           ye       zeros(n.cc)     zeros(n.cc);
            zeros(n.cc)                ye     zeros(n.cc);
            zeros(n.cc)       zeros(n.cc)             ye];
  
end

% Formulação geral
Yt = Ye;

%% INCLUSÃO DO EFEITO DA ADMITÂNCIA DO SOLO: ELEMENTOS MÚTUOS DA MATRIZ Y
% Yg=zeros(n.cc,n.cc);
FlagPgAprox=1;
if FlagZg==1; FlagPgAprox=0; end

[Pg] = SCC_Admitancia(n,p,s,rho.g,er.g,FlagPgAprox);

% if n.f==1 
%     Yg(end) = yg;
% end
% 
% % Formulação geral
% Yt = (Ye\eye(n.c)+Yg\eye(n.c))\eye(n.c);

end
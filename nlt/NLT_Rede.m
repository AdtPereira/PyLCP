function [V] = NLT_Rede(s,Y,GAM,TI,Lx,R,L,T,ksi,n,FlagRd)

%% Declaração das variáveis
YS=zeros(n.c,n.c); YL=zeros(n.c,n.c); Nbc=n.b*n.c; 
V=zeros(Nbc,n.N); I=zeros(Nbc,1); VSp=1; Vsrc=1; Xopt=0;

%% FONTE DE TENSÃO
if Vsrc==1    
    tdelay=ksi*T; VS=(VSp./s).*exp(-tdelay*s);
elseif Vsrc==2
    vt=VSp.*ones(1,n.N).*exp(-c*ts); vt(1)=0; vt(round(twidth/dt):n.N)=0; 
    VS=dt*fft(vt);
elseif Vsrc==3
    wt = 2*pi*60; VS = VSp*s./(s.^2 + wt^2);
end

%% TOPOLOGIA DE REDE 
if Xopt == 0    
    ZS = R.S + s*L.S;       % Impedância da fonte no terminal emissor    
    ZL = R.L + s*L.L;   	% Impedância da carga no terminal receptor
    ZG = R.G + s*L.G;   	% Impedância da carga no terminal receptor
else
    fss = Xopt;
    ZS = R.S + s*(L.S/(2*pi*fss)); % Impedância da fonte terminal emissor
    ZL = R.L + s*(L.L/(2*pi*fss)); % Impedância da carga terminal receptor
    ZG = R.G + s*(L.G/(2*pi*fss)); % Impedância da carga terminal receptor
end

%% FORMULAÇÃO MATRICIAL
for k=1:n.Nr            
    %% MATRIZES DE REDE
                        % Terminais Emissor e Receptor a vazio
                        % MODELO M11SCC1C
                        % MODELO M21SCC1C 
                        % MODELO M31SCC2C
                        % MODELO M41SCC2C
                        % MODELO M91SCC3C
                        YS(1,1)=1/ZS(k);

    if FlagRd==3        % Terminais Emissor aterrados
                        % MODELO M112SCC2C  (MORCHED ET AL.. 1999)
                        YS(2,2)=1/ZG(k); 
                        YS(3,3)=1/ZG(k); 
                        YS(4,4)=1/ZG(k); 

    elseif FlagRd==5    % Terminal Emissor do segundo condutor Aterrado
                        % MODELO M52SCC1C
                        % MODELO M73SCC1C   (DUARTE, 2022)
                        % MODELO M101SCC2C  (MORCHED ET AL., 1999)
                        % MODELO M121SCC2C  (MORCHED ET AL., 1999)
                        YS(2,2)=1/ZG(k);      

    elseif FlagRd==6    % Terminal Emissor do segundo condutor Aterrado
                        % Terminal Emissor do terceiro condutor Aterrado
                        % MODELO M63SCC1C
                        YS(2,2)=1/ZG(k); 
                        YS(3,3)=1/ZG(k); 
    end

    %% REPRESENTAÇÃO MATRICIAL DA LINHA DE TRANSMISSÃO
    Ti=TI(:,:,k); Gam=GAM(:,:,k); 
    
    % Impedância característica no domínio das fases a partir 
    % da segunda representação de quadripolos obtida no domínio modal
    % Ver Item 7.8 (DE CONTI, 2023, p. 63)
    Yc=Ti/Gam/Ti*Y(:,:,k);
    
    % Definição Coth Matricial
    cothm = Ti/(expm(Gam*Lx) - expm(-Gam*Lx))*...
               (expm(Gam*Lx) + expm(-Gam*Lx))/Ti;
    
    % Definição csch Matricial
    cschm = 2*Ti/(expm(Gam*Lx) - expm(-Gam*Lx))/Ti;
    
    % Representação Matricial da LT
    Ykk=cothm*Yc; Ykm=-cschm*Yc;     
        
    %% CORRENTES INJETADAS NAS BARRAS
    I(1)=VS(k)/ZS(k);

    %% FORMULAÇÃO DO PROBLEMA V = ZI (Ax = b)
    V(:,k) = [YS+Ykk Ykm;Ykm Ykk+YL]\I;
end
end
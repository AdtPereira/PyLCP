function [Z,Y,Zc,Yc,Gamma] = ohtl_Parametros_QuaseTEM(jw,h,ro,ri,sigc,rhog,erg,FLAG,APLOG)
addpath 'C:\Users\adilt\OneDrive\01 ACADEMIA\06 MODELOS\LCP\basic'

%% 1. CONSTANTES FÍSICAS
if FLAG == 1; erg=1; end    % Aproximação de CARSON
eo = 8.854E-12;             % Permissividade absoluta do vácuo
uo = 4*pi*1E-7;             % Permeabilidade magnética vácuo 
eg = erg*eo;                % Permissividade do solo
sigg = 1/rhog;              % Condutividade elétrica do solo [S/m]

% Condutividade elétrica do condutor [S/m]
% sigc = 1./(rcc.*pi.*(ro.^2 - ri.^2));  
rcc = 1./(sigc.*pi.*(ro.^2 - ri.^2));

%% 2. DECLARAÇÃO DAS VARIÁVEIS
Nc = length(h); Ns = length(jw); Z = zeros(Nc,Nc,Ns); Y = zeros(Nc,Nc,Ns);
Zc = zeros(Nc,Nc,Ns); Yc = zeros(Nc,Nc,Ns); Gamma = zeros(Nc,Nc,Ns);
Zi = zeros(Nc,Nc,Ns); S1 = zeros(Nc,Nc); S2 = zeros(Nc,Nc); 
T = zeros(Nc,Nc); COND = zeros(Nc,3);

%% 3. MATRIZ AUXILIAR M
for k=1:Nc;    COND(k,:) = [0 h(k) ro(k)]; end 
[MD, Md, Mx] = FunMatrizDistancias(COND); M = log(MD./Md);

for k=1:Ns
    %% 4. IMPEDÂNCIA INTERNA, Zi
    zihf = 1./(2*pi.*ro).*sqrt(uo*jw(k)./sigc); 
    zi = sqrt(rcc.^2 + zihf.^2);
    for i=1:Nc
        Zi(i,i,k) = zi(i); 
    end

    %%  5. Matriz de Admitância com retorno pelo solo, Yg

    y1 = jw(k)*sqrt(uo*eo);    	% Constante de propagação intrínsica do 
                                % meio 1 (ar)
                                % y1 [Pg14-Cp4]
                                
    y2 = sqrt(jw(k)*uo*(sigg + jw(k)*eg));
                                % Constante de propagação intrínsica do 
                                % meio 2 (solo)
                                % y1 [Pg14-Cp4]
    n = y2/y1;                  % Coef. de refração do solo [Pg25-Cp4]
    
    if APLOG
        %% 5.1 Constante de propagação na direção transversal [Pg9-Cp4]  
        ni = sqrt(y2^2 - y1^2);   
        for i=1:Nc
            for j=1:Nc 
                Num = n^2+1;
                Den = ni*sqrt((h(i)+h(j))^2 + (Mx(i,j))^2);
                
                S1(i,j) = log(1+2/Den);  
                S2(i,j) = 2/Num*log(1+Num/Den);
                T(i,j) = 2*log(2) + 2*n^2/Num*log((1+Num/Den)/(1+2*Num/Den));
            end
       end
    else
       for i=1:Nc
           for j=1:Nc            
                F1 = @(x) exp(-x*(h(i)+h(j)))/...
                        (x + sqrt(x^2 + y2^2 - y1^2)) * cos(Mx(i,j)*x);    
                S1(i,j) = 2 * integral(F1,0,Inf,'ArrayValued',true);

                F2 = @(x) exp(-x*(h(i)+h(j)))/...
                        (n^2*x + sqrt(x^2 + y2^2 - y1^2)) * cos(Mx(i,j)*x);
                S2(i,j) = 2 * integral(F2,0,Inf,'ArrayValued',true);

                F3 = @(x) sqrt(x^2 + y2^2 - y1^2)/x * ...
                    (exp(-0.5*x*(h(i)+h(j))) - exp(-1.0*x*(h(i)+h(j))))/...
                        (n^2*x + sqrt(x^2 + y2^2 - y1^2)) * cos(Mx(i,j)*x);
                T(i,j) = 2 * integral(F3,0,Inf,'ArrayValued',true);   
            end
       end
    end

    if FLAG == 1        
        %% 5.2 APROXIMAÇÃO DE CARSON [Pg18-Cp4]
        % Adotando-se que a tensão V(x) corresponda ao potencial elétrico
        % na superfície do condutor em relação ao infinito: [Pg19-Cp4]   
        Z(:,:,k) = Zi(:,:,k) + jw(k)*uo/2/pi*(M + S1);
        Y(:,:,k) = jw(k)*2*pi*eo*inv(M);  
        
     elseif FLAG == 2   
        %% 5.3. APROXIMAÇÃO DE SUNDE [Pg19-Cp4]
        % Na expressão proposta por Nakagawa (1981), Se o termo (𝜀𝑟𝑔 − 1)
        % for substituído  simplesmente por  𝜀𝑟𝑔, obtém-se a equação ,
        % de Sunde.   
        Z(:,:,k) = Zi(:,:,k) + jw(k)*uo/2/pi*(M + S1);
        Y(:,:,k) = jw(k)*2*pi*eo*inv(M); 
        
     elseif FLAG == 3
        %% 5.4. APROXIMAÇÃO DE NAKAGAWA [Pg19-Cp4]
        % Adotando-se que a tensão V(x) corresponda ao potencial elétrico
        % na superfície do condutor em relação ao infinito: [Pg19-Cp4] 

        % Equação integral
        for i=1:Nc
           for j=1:Nc
                F1 = @(x) exp(-x*(h(i)+h(j)))/...
                        (x + sqrt(x^2 + jw(k)*uo*(sigg+jw(k)*(erg-1)*eo))) * cos(Mx(i,j)*x);    
                S1(i,j) = 2 * integral(F1,0,Inf,'ArrayValued',true);
        
                F2 = @(x) exp(-x*(h(i)+h(j)))/...
                        (n^2*x + sqrt(x^2 + jw(k)*uo*(sigg+jw(k)*(erg-1)*eo))) * cos(Mx(i,j)*x);
                S2(i,j) = 2 * integral(F2,0,Inf,'ArrayValued',true);
           end
        end

        Z(:,:,k) = Zi(:,:,k) + jw(k)*uo/2/pi*(M + S1);
        Y(:,:,k) = jw(k)*2*pi*eo*inv(M + S2);  
        
    elseif FLAG == 4   
        %% 5.5. APROXIMAÇÃO QUASE-TEM (REP. INTEGRAL) [Pg24-Cp4]
        % Supõe-se, por aproximação, que a constante de propagação da 
        % corrente seja igual à constante de propagação no ar.  
        
        % SOLUÇÃO PARA A DEFINIÇÃO DE TENSÃO (A)
        % Supõe-se que r << h, escreve-se uma solução quase-TEM  
        % considerando que a tensão V(x) corresponde à integral do campo  
        % elétrico vertical da superfície do solo à altura do condutor  
        % incluindo a componente z do potencial vetor magnético. [Pg17-Cp4]    
        Z(:,:,k) = Zi(:,:,k) + jw(k)*uo/2/pi*(M + S1 - (T + S2));
        Y(:,:,k) = jw(k)*2*pi*eo*inv(M - T);
    end
    
    %% 6. RESULTADOS PARCIAIS DA APROXIMAÇÃO QUASE-TEM
    Zc(:,:,k) = Y(:,:,k)\sqrtm(Y(:,:,k)*Z(:,:,k));
    Yc(:,:,k) = Z(:,:,k)\sqrtm(Z(:,:,k)*Y(:,:,k));
    Gamma(:,:,k) = sqrtm(Z(:,:,k)*Y(:,:,k)); 
end
end
function [Z,Ye,Zc,Yc,TAUv] = ohtl_Parametros_Carson(h,ro,ri,rcc,rhog,erg,Ge,s)
%CalculaParametrosLog Summary of this function goes here
%   Detailed explanation goes here

%% Constantes Físicas
  eo = 8.854*1E-12; % Permissividade absoluta do vácuo
  uo = 4*pi*1E-7;   % Permeabilidade magnética vácuo
sigg = 1/rhog;      % Condutividade elétrica do solo [S/m]

%% Declaração das variáveis auxiliares
Nc = length(h);   
Ns = length(s); 
COND = zeros(Nc,3);

%% Cálculo das Matrizes Auxiliares de Distânci
for k=1:Nc
    COND(k,:) = [0 h(k) ro(k)];
end
[MD, Md, Mx] = FunMatrizDistancias(COND);

%% 1. Fórmula aproximada da Resistência interna em corrente contínua de 
% condutores cilíndricos tubulares
% Ref.: Notas de aulas. Capítulo 3, Item 2.2.B, página 6. (DE CONTI, 2023)
% 1. O Cabo de alumínio CAA será modelado como um condutor cilíndrico 
% tubular;
%
% 2. Dado que a condutividade do condutor varia com a frequência, seu valor 
% será calculado a partir do valor da resistência em corrente contínua do
% condutor na temperatura desejada.

sigc = 1./(rcc.*pi.*(ro.^2 - ri.^2));
zihf = 1./(2*pi.*ro).*sqrt(uo*s./sigc);
zi = sqrt(rcc.^2 + zihf.^2);

Zi = zeros(Nc,Nc,Ns);  
for k=1:Ns
    for l=1:Nc
        Zi(l,l,k) = zi(l,k);
    end
end

%% 2. Matriz de impedâncias externas 
% Ref.: Notas de aulas. Capítulo 3, Item 2.3, página 7. (DE CONTI, 2023)
% 1. A matriz de indutância externa, Le, é invariante com a frequência.

Le = uo/2/pi*log(MD./Md);
Ze = zeros(Nc,Nc,Ns);

for k=1:Ns
    Ze(:,:,k)=s(k)*Le; 
end

%% 3. Matriz de impedâncias de retorno no solo 
% Ref.: Notas de aulas. Capítulo 3, Item 2.4.D, pág. 11. (DE CONTI, 2023)
% 1. A matriz de impedâncias com retorno no solo é calculada pelas 
% aproximações logarítmicas propospostas por Deri et al. (1981).

yg = sqrt(s*uo.*(sigg+s*erg*eo));
 p = 1./yg;
zg = zeros(Nc,Nc,Ns);
Zg = zeros(Nc,Nc,Ns);

for k=1:Ns
    for l=1:Nc
        for m=1:Nc
            if l==m
                zg(l,m,k) = log((h(l)+p(k))/h(l));
            else
                zg(l,m,k) = log(sqrt((h(l)+h(m)+2*p(k))^2 + Mx(l,m)^2)/...
                               (sqrt((h(l)+h(m)+0*p(k))^2 + Mx(l,m)^2)));     
            end
        end
    end
    Zg(:,:,k) = s(k)*uo/2/pi*zg(:,:,k);
end

%% 4. Matriz de Impedâncias total
Z = Zi + Ze + Zg;

%% 5. Matriz de Admitância Externa, Ye, e Impedância Característica, Zc
% Ref.: Notas de aulas. Capítulo 3, Item 3.1, pág. 13. (DE CONTI, 2023)
% 1. A matriz de capacitâncias externas, Ce, é invariante com a frequência.
Ce = uo*eo*inv(Le);         

Zc = zeros(Nc,Nc,Ns);
Yc = zeros(Nc,Nc,Ns);       
Ye = zeros(Nc,Nc,Ns); 
TAUv = zeros(Nc,Nc,Ns);       

for k=1:Ns
    Ye(:,:,k) = eye(Nc)*Ge + s(k)*Ce;
    Zc(:,:,k) = Ye(:,:,k)\sqrtm(Ye(:,:,k)*Z(:,:,k));
    Yc(:,:,k) = Z(:,:,k)\sqrtm(Z(:,:,k)*Ye(:,:,k));
    TAUv(:,:,k) = sqrtm(Z(:,:,k)*Ye(:,:,k));
end
end
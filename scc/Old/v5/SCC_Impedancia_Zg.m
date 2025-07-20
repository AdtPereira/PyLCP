function [Zg] = SCC_Impedancia_Zg(n,p,s,rhog,erg,Flag)
% UNTITLED Summary of this function goes here
% Detailed explanation goes here

%% CONSTANTES FÍSICAS E DADOS DE ENTRADA 
eo = 8.854E-12;                    	% Permissividade absoluta do vácuo
uo = pi*4E-7;                     	% Permeabilidade magnética vácuo

sg = 1/rhog;                      	% Condutividade do solo
jw = s;                  	        % Frequência angular

hij = p(:,2);                   	% Profundidade dos condutores
Zg = zeros(n.f,n.f);                % Matriz Imped. Retorno pelo solo                                     

[Dij,dij,xij] = SCC_Distancias(p);  % Matrizes distâncias relativas

%% CONSTANTES DE PROPAGAÇÃO
                                    % FLAG 1 A. MAGALHAES/XUE (2015)
yo = jw*sqrt(uo*eo);          	    % Constante de Propagação do ar
yg = sqrt(jw*uo*(sg + jw*erg*eo));	% Constante de Propagação do Solo
                                    
if Flag == 2; yo = 0; end           % FLAG 2 B. FORMULAÇÃO DE SUNDE 

if Flag == 3 || Flag == 4 || Flag == 6 || Flag == 8  
                                    % APROXIMAÇÕES EM BAIXAS FREQUÊNCIAS
                                    % FLAG 3 C.     POLLACZEK (1931)    
                                    % FLAG 4 D.     AMETANI (1931)
                                    % FLAG 6 A.2    SAAD ET AL. (1996)
                                    % FLAG 8 A.4    WEDEPOHL/WILCOX (1973)    
    yo = 0; yg = sqrt(jw*uo*sg);         	
end

%% A.4 APROXIMAÇÃO DE WEDEPOHL E WILCOX (1973)
if Flag == 8    
    % Uma  aproximação para a expressão de Pollaczek foi proposta por 
    % Wedepohl e Wilcox (1973) para o caso em que |𝑚𝑑| < 0,25:    
    m = yg; y = 0.577215665;
    for i=1:n.f
        for j=1:n.f
            Zg(i,j) = jw*uo/2/pi*...
                (-log(y*m*dij(i,j)/2) + 0.5 - 2/3*m*(hij(i)+hij(j)));    
        end
    end
end

%% A.3 APROXIMAÇÃO DE THEETHAYI ET AL. (2007)
if Flag == 7 
    % Uma aproximação logarítmica para as fórmulas de Sunde foi  
    % proposta por Theethayi et al. (2007), sem prova matemática, para 
    % o caso  especial em que os cabos estão na mesma profundidade,
    % ℎ𝑖 = ℎ𝑗:
    for i=1:n.f
        for j=1:n.f
            Zg(i,j) = jw*uo/2/pi*(log((1 + yg*xij(i,j))/yg/xij(i,j)) ...
                  + exp(-(hij(i) + hij(j))*yg)*(2/(4 + (yg*xij(i,j))^2)));    
        end
    end
end

%% A.1 APROXIMAÇÃO DE DE CONTI ET AL. (2023)
if Flag == 5 || Flag == 6    
    % FLAG == 5 - A expressão foi proposta por De Conti et al. (2023) como 
    % aproximação para a equação integral de Magalhães/Xue.    
    
    % FLAG == 6 - A.2 APROXIMAÇÃO DE SAAD ET AL. (1996)
    % Se na aproximação de De Conti et al. (2023) for considerado 
    % 𝛾0 = 0,  obtém-se uma aproximação para a equação de Sunde.
    % No caso particular em que 𝛾𝑔 = sqrt(𝑗𝜔𝜇o𝜎𝑔), essa expressão se 
    % torna igual à expressão proposta por Saad et al. (1996) como  
    % aproximação para a fórmula de Pollaczek.  
    for i=1:n.f
        for j=1:n.f
            Zg(i,j) = jw*uo/2/pi*(besselk(0,yg*dij(i,j)) + ...                      
                        (yg-yo)/(yo+yg)*exp(-(hij(i) + hij(j))*yg)* ...
                            (2/(4 + (yg*xij(i,j))^2)));    
        end
    end
end 
 
%% A. EQUAÇÕES INTEGRAIS
 
if Flag == 4  
    % D. FORMULAÇÃO DE AMETANI (1931)
    % Consideração 1. 𝜎𝑔 ≫ 𝜔𝜀𝑟𝑔𝜀o. [Pg24-Cp5]
    % Consideração 2. Aproximação da integral de Pollaczek pela  
    %                 integral de Carson, o que corresponde a fazer com 
    %                 que 𝛾𝑔 = 0 no coeficiente do termo exponencial 
    %                 presente no numerador dessa integral. [Pg25-Cp5]    
    for i=1:n.f
        for j=1:n.f
            F1 = @(x) cos(xij(i,j)*x)*exp(-(hij(i) + hij(j)) * x)/...
                         (x + sqrt(x^2 + yg^2));
            
            Zg(i,j) = jw*uo/2/pi*...
                (besselk(0,yg*dij(i,j)) - besselk(0,yg*Dij(i,j))+... 
                        2*integral(F1,0,Inf,'ArrayValued',true));         
        end
    end
end

if Flag == 1 || Flag == 2 || Flag == 3  
    % A. FORMULAÇÃO DE MAGALHAES/XUE ET AL.(2018)
    % A expressão mais rigorosa para o cálculo impedância de retorno do 
    % solo de cabos enterrados foi proposta inicialmente por Magalhães e 
    % Lima (2015) e posteriormente por Xue et al. (2018). Sua representação
    % mais compacta tem a forma (De Conti et al., 2023):
    
    % B. FORMULAÇÃO DE SUNDE
    % A impedância de retorno do solo de um cabo enterrado proposta por 
    % Sunde é um caso particular da equação de Magalhães/Xue se 𝛾0 = 0.
    
    % C. FORMULAÇÃO DE POLLACZEK (1931) 
    % No caso em que 𝜎𝑔 ≫ 𝜔.𝜀𝑟𝑔.𝜀0, a expressão integral proposta por Sunde 
    % se reduz às expressões clássicas proposta por Pollaczek (1931).
    % A fórmula de Pollaczek é uma aproximação válida somente em baixas 
    % frequências e solos de baixa resistividade, sendo análoga à expressão
    % desenvolvida por Carson para linhas aéreas.    
    for i=1:n.f
        for j=1:n.f
            F2 = @(x) cos(xij(i,j)*x)*exp(-(hij(i) + hij(j))...
                * sqrt(x^2 + yg^2))/(sqrt(x^2 + yo^2) + sqrt(x^2 + yg^2));
                  
            Zg(i,j) = jw*uo/2/pi*(besselk(0,yg*dij(i,j))...
                - besselk(0,yg*Dij(i,j))+... 
                        2*integral(F2,0,Inf,'ArrayValued',true));
        end
    end 
end
end    
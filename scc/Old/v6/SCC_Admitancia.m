function [Pg,Yg] = SCC_Admitancia(n,p,s,rhog,erg,FlagPgAprox)
% UNTITLED Summary of this function goes here
% Detailed explanation goes here

%% CONSTANTES FÍSICAS E DADOS DE ENTRADA 
eo = 8.854E-12;                    	% Permissividade absoluta do vácuo
uo = pi*4E-7;                     	% Permeabilidade magnética vácuo
eg = erg*eo;                     	% Permissividade do solo
sg = 1/rhog;                      	% Condutividade do solo
jw = s;                  	        % Frequência angular
hij = p(:,2);                   	% Profundidade dos condutores
[Dij,dij,xij] = SCC_Distancias(p);  % Matrizes distâncias relativas

%% CONSTANTES DE PROPAGAÇÃO
yo = jw*sqrt(uo*eo);              	% Constante de Propagação do ar
yg = sqrt(jw*uo*(sg + jw*eg));      % Constante de Propagação do Solo

%% ADMITÂNCIA ASSOCIADA AO RETORNO PELO SOLO, Yg
Pg=zeros(n.f,n.f); 

%%  Formulação de Xue para um solo dispersivo. [Pg31-Cp5]
for i=1:n.f
    for j=1:n.f
        if FlagPgAprox
            %%  Expressão aproximada De Conti et al. (2023) [Pg32-Cp5]
            Pg(i,j) = jw/2/pi/(sg + jw*eg)*...
                (besselk(0,yg*dij(i,j)) + (yg^2-yo^2)/(yg^2+yo^2)* ...
                    besselk(0,yg*Dij(i,j)));          

        else
            F = @(x) sqrt(x^2+yo^2)/sqrt(x^2+yg^2)* ...
                    exp(-(hij(i)+hij(j))*sqrt(x^2+yg^2))*cos(xij(i,j)*x)/...
                        (sqrt(x^2+yo^2) + yo^2/yg^2*sqrt(x^2+yg^2));
    
            Pg(i,j) = jw/2/pi/(sg + jw*eg)*...
                        (besselk(0,yg*dij(i,j)) - besselk(0,yg*Dij(i,j)) + ... 
                        2*integral(F,0,Inf,'ArrayValued',true));            
        end
    end
end

% Admitância de retorno de solo
Yg = jw*Pg\eye(n.f);

end
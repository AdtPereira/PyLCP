function [Ye] = SCC_Admitancia_BaixasFreq(scc,s,er,n)

% UNTITLED Summary of this function goes here
% Detailed explanation goes here

% Desprezar a admitância do solo significa supor que o solo se comporte 
% como um condutor elétrico perfeito. Isso restringe a aplicabilidade da 
% formulação matricial acima a situações que envolvam solos com baixa 
% resistividade. Em solos de alta resistividade ou no estudo de fenômenos 
% de alta frequência, a matriz acima deve ser modificada de forma a incluir
% a admitância do solo.

% A rotina de cálculo de parâmetros de cabos disponível no ATP despreza a 
% admitância do solo. Com isso, supõe que os cabos estejam imersos em um 
% meio condutor ideal.

% Desprezar-se 𝒀𝒈′ pode ser considerado uma boa aproximação somente em 
% estudos de fenômenos de baixas frequências em solos de baixa 
% resistividade.

%% CONSTANTES FÍSICAS
eo=8.854E-12;           % Permissividade absoluta do espaço livre
jweo=s*eo;   	        % Frequência angular

%% ADMITÂNCIA TRANSVERSAL EM BAIXAS FREQUÊNCIAS
Y1 = jweo*2*pi*er.ip/log(scc.b/scc.a); % Admitância da isolação Primária
Y2 = jweo*2*pi*er.is/log(scc.d/scc.c); % Admitância da isolação Secundária
Y3 = jweo*2*pi*er.ie/log(scc.f/scc.e); % Admitância da isolação Externa 

%% MATRIZ ADMITÂNCIA EXTERNA
if n.cc==1
    Ye = Y1;
elseif n.cc==2
    Ye = [Y1 -Y1;-Y1 Y1+Y3];
elseif n.cc==3
    Ye = [Y1 -Y1 0;-Y1 Y1+Y2 -Y2;0 -Y2 Y2+Y3];
end 
end
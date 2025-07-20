function [Ti,Lmbd,NrIte] = Mdal_EigNR(Z,Y,Ti,Lmbd,s,EpsTol,IterMax)
% Rotina para garantia da variação contínua dos autovetores associados 
% ao produto matricial YZ, no domínio da frequência, usando a técnica de 
% Newton-Raphson (NR).
%
% [1]   WEDEPOHL, L. M, H. V. Nguyen, G. D. Irwin, "Frequency-dependent 
%       transformation matrices for untransposed transmission lines 
%       using Newton-Raphson method", IEEE Trans. Power Systems, 
%       vol. 11, no. 3, 1538-1546, Aug, 1996.
%
% [2]   ZLATUNIĆ, I.; VODOPIJA, S.; GOIĆ, R., "Mode switching in modal 
%       domain models of overhead lines and underground cables", 
%       International Conference on Power Systems Transients (IPST), 2015,
%       ISSN: 2434-9739. https://www.ipstconf.org/Proc_IPST2015.php
%
% [3]   CHRYSOCHOS, Andreas I.; PAPADOPOULOS, Theofilos A.; PAPAGIANNIS, 
%       Grigoris K., "Robust Calculation of Frequency-Dependent 
%       Transmission-Line Transformation Matrices Using the 
%       Levenberg–Marquardt Method", IEEE Trans. Power Delivery, 
%       vol. 29, no. 4, 162-1629, Aug, 2014.
%
%                                           CONTAGEM, 26 DE JANEIRO DE 2024

% EIGENPROBLEM SOLVING WITH NEWTON-RAPNSON (WEDEPOHL, 1996)
% S         Produto Y*Z
% Tik       Coluna k da Matriz de Transformação de corrente, Ti: Autovetor 
%           associado ao condutor k
% lkk       Autovalor associado ao condutor k
% x         Vetor coluna incognita: Tik e o autovalor lkk
% F         Função polinomial
% J         Matriz Jacobiana 

% Variáveis Auxiliares
Nc=size(Z,1); U=eye(Nc); NrIte=0; 

% Normalização Matricial (Item 7.2 de [1])
eo=8.854E-12; uo=pi*4E-7; w=imag(s); Norm=-w^2*eo*uo;
S=Y*Z/Norm-U;  Lmbd=Lmbd/Norm-U;       

for k=1:Nc 
    % Inicializa as variáveis
    Stop=0;
    Tik=Ti(:,k); lkk=Lmbd(k,k); x=[Tik; lkk];                               
    F=[(S-lkk*U)*Tik; sum(Tik.^2)-1];           
    J=[S-lkk*U; 2*Tik.']; J(:,Nc+1)=[-Tik;0];   

    % Inicializa o processo interativo de NR
    while ~Stop 
        % Resolve o problema de NR
        x=x-J\F;

        % Atualiza os autovetores e autovalores
        Tik=x(1:Nc); lkk=x(Nc+1);   
        
        % Atualiza a Função F
        F=[(S-lkk*U)*Tik; sum(Tik.^2)-1]; 
        
        % Critério de Convergência - Ver Item 2.B de [3]
        if sum(abs(F))<EpsTol; Stop=1; end         

        % Verifica as Interações
        NrIte=NrIte+1; if NrIte==IterMax; Stop=1; end              
            
        % Atualiza Matriz Jacobiana    
        J=[S-lkk*U; 2*Tik.']; J(:,Nc+1)=[-Tik;0]; 
    end
    Ti(:,k)=Tik;            % Atualiza Matriz de Transformação de corrente 
    Lmbd(k,k)=lkk;
    Lmbd(k,k)=Norm*(1+lkk); % Atualiza autovalor associado ao condutor k
end
end
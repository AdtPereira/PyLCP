function [R] = ReordernaMatriz(n,R)

Op1=1:n.cc:n.c;                     % Define operações

                                    % Realiza as trocas dos elementos 
                                    % próprios dos condutores n.cc=1
for k=1:n.f
    R = TrocarPos(R,Op1(k),k,1);    % Operações nas Linhas de R
    R = TrocarPos(R,Op1(k),k,2);    % Operações nas Colunas de R
end
                                    % Realiza as trocas dos elementos 
                                    % próprios das blindagens n.cc=2
    R = TrocarPos(R,n.f+1,n.f+2,1); % Operações nas Linhas de R
    R = TrocarPos(R,n.f+1,n.f+2,2); % Operações nas Colunas de R
end
            
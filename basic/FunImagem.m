    function X1 = FunImagem (X)
    % Imagem Calcula a imagem de um vetor bidimensional com relação ao 
    % eixo horizontal.
    % A função Imagem calcula a imagem do vetor X(Xx, Xy) com relação ao 
    % eixo horizontal y=0 no sistema de coordenadas cartesiano.
    %
    % Calling sequence:
    % Xi = Imagem (X)
    %
    % Definição de variáveis
    % X             -- Posição original.
    % X1            -- Imagem horizontal da posição original.
    %
    % Memória de revisão:
    % Data            Programador                  Descrisão da mudança
    % =====           ==============               ====================
    % 24-03-2018      Adilton Junio L. Pereira     Implementação TCC-II
    %
    % Cálculo efetivo da distancia
    X1 = [X(1), -X(2)];
    end
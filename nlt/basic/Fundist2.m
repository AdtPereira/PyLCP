    function distance = Fundist2 (X1, X2)
    %DIST2 calculate the distance between two points with reference.
    %Function DIST2 calculates the distance between two points X1 (x1,y1) 
    %and X2 (x2,y2) in a Cartesian coordinate system.
    %
    %   Calling sequence:
    %   distance = dist2 (X1, X2)

    %Definição de variáveis
    %x1             - Abscissa do ponto X1. Primeira Posição do vetor X1.
    %y1             - Ordenada do ponto X1. Segunda  Posição do vetor X1.
    %x2             - Abscissa do ponto X2. Primeira Posição do vetor X2.
    %y2             - Ordenada do ponto X2. Segunda  Posição do vetor X2.
    %ox             - Abscissa do sistema de origem. 
    %                 Primeira Posição do vetor O.
    %distance       - Distância entre os dois pontos.

    % Memória de revisão:
    %   Data            Programador                    Descrisão da mudança
    %   =====           ==============                 ====================
    %   02-01-2007      S. J. Chapman                  Código original
    %   17-08-2013      Adilton Junio L. Pereira       Implem. para estudo
    %   24-03-2018      Adilton Junio L. Pereira       Implem. para TCC-II
    %   02-12-2020      Adilton Junio L. Pereira       Rev. Dissertação

    %Cálculo efetivo da distancia
    distance = power((X2(1)-X1(1))^2+(X2(2)-X1(2))^2,0.5);
    end
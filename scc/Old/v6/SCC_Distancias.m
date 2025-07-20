function [Dij, dij, xij] = SCC_Distancias(M)
    %Calcula as matrizes de posição relativa entre os
    %cabos coaxiais subterrâneos do tipo singe core
    %
    %   DADOS DE ENTRADA
    %
    %   M Matriz posição com par cartesiano dos condutores e raio 
    %     externo do cabo coaxial.
    %    
    %   M = [ xi yi ri_ext ] 
    %       [ .. .. .. ... ]
    %       [ xj yj rj_ext ]
    %    
    %   onde 
    %         xi = Cota x do cond. i, em relação à referência
    %         yi = hi, profundidade positiva do cond. i, em relação ao solo
    %     ri_ext = Raio externo do cond. i
    %
    %   DADOS DE SAÍDA
    %   MD é a distância entre o condutor i a imagem do condutor j
    %   MD = sqrt((hi+hj)^2 + xij^2)
    %   MD = [ Dii ... Dij ] 
    %        [ ... ... ... ]
    %        [ Dji ... Djj ]
    %                        
    %   Md = [ dii ... dij ] 
    %        [ ... ... ... ]
    %        [ dji ... djj ]
    %                        
    %   Mx = [ dii ... xij ] 
    %        [ ... ... ... ]
    %        [ xji ... djj ]
    %                        
    %   onde 
    %        dii = ri
    %        djj = rj
    %         yi = hi   Altura do cond.-i, em relação ao solo
    %         ri        Raio do cond.-i 
    %
    
    nc=size(M,1); Dij=zeros(nc,nc); dij=zeros(nc,nc); xij=zeros(nc,nc);
                                              
    for i = 1:nc
        for j = 1:nc
            if i == j 
                xij(i,j) = M(i,3);
            else
                xij(i,j) = abs(M(j,1) - M(i,1));
            end
            dij(i,j) = sqrt((M(i,2) - M(j,2))^2 + xij(i,j)^2);
            Dij(i,j) = sqrt((M(i,2) + M(j,2))^2 + xij(i,j)^2);
        end           
    end
end
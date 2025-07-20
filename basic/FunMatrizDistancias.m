function [MD, Md, Mx] = FunMatrizDistancias(COND)
    %CalcularDistancias Calcula as matrizes de posição relativa entre os
    %cabos coaxiais subterrâneos do tipo singe core
    %
    %   DADOS DE ENTRADA
    %
    %   COND -  Matriz posição com par cartesiano dos condutores e raio 
    %           externo do cabo coaxial.
    %    
    %   COND = [ xi yi ri_ext ] 
    %          [ .. .. .. ... ]
    %          [ xj yj rj_ext ]
    %    
    %   onde 
    %         xi = Cota x do cond. i, em relação à referência
    %         yi = hi, profundidade positiva do cond. i, em relação ao solo
    %     ri_ext = Raio externo do cond. i
    %
    %   DADOS DE SAÍDA
    %                            
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
    %         yi = hi   Altura do cond. -i, em relação ao solo
    %         ri        Raio do cond. -i 
    %

    [NC, ~] = size(COND);
         MD = zeros(NC,NC); 
         Md = zeros(NC,NC);
         Mx = zeros(NC,NC);

    for k = 1:NC
        for l = 1:NC
           MD(k,l) = Fundist2(COND(k,:), FunImagem(COND(l,:)));   
           Md(k,l) = Fundist2(COND(k,:), COND(l,:));
           Md(k,k) = COND(k,3);
           Mx(k,l) = abs(COND(l,1) - COND(k,1));
           Mx(k,k) = COND(k,3);
        end
    end
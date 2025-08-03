function [M] = TrocarPos(M,a,b,Index)

if Index == 1 % Altera Linhas
    La = M(a,:); Lb = M(b,:); M(a,:) = Lb; M(b,:) = La;
elseif Index == 2 % Altera Colunas
    Ca = M(:,a); Cb = M(:,b); M(:,a) = Cb; M(:,b) = Ca;
else
    fprintf('\n\n\n Operação não Permitida. Matriz inalterada! \n\n');

end
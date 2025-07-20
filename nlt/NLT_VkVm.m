function [Vk,Vm] = NLT_VkVm(Z,Y,Lx,VS,RS,RL,s)
%Dados de Entrada
%
% Z: Impedância longitudinal, ohm/m
% Y: Admitância transversal, S/m
% Lx: Comprimento da linha, m
% VS: Tensão aplicada no domínio da frequência, V
% RS: Resistência interna da fonte
% RL: Resistência da carga
% s: Vetor de frequências complexas (s+i*w)

% Variáveis de Saída
%
% Vk: Tensão no terminal emissor, V
% Vm: Tensao no terminal receptor, V

% Parâmetros da linha
gama=sqrt(Z.*Y);
Yc=sqrt(Y./Z);

N=length(s);
V(2,N)=0;

for n=1:N
    Ykk=Yc(n)*coth(gama(n)*Lx);
    Ymm=Ykk;
    Ykm=-Yc(n)*csch(gama(n)*Lx);
    %Ykm=-Yc(n)*2/(exp(gama(n)*Lx)-exp(-gama(n)*Lx));
    Ymk=Ykm;
    Ybus=[Ykk+1/RS  Ykm;
          Ymk       Ymm+1/RL];
      V(:,n)=Ybus\[VS(n)/RS; 0];
end

Vk=V(1,:);
Vm=V(2,:);
end


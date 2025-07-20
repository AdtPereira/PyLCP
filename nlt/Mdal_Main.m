function [GAM,TI,TIEIG,Zcm,Ups,NrIte] = Mdal_Main(Z,Y,s,n,IterMax,EpsTol) 

% Declara variáveis iniciais nulas
TI=zeros(n.c,n.c,n.Nr); GAM=zeros(n.c,n.c,n.Nr); TIEIG=zeros(n.c,n.c,n.Nr);
NRIte=zeros(1,n.Nr); Zcm=zeros(n.c,n.c,n.Nr); Ups=zeros(n.c,n.c,n.Nr);

% Define Valor Inicial 
[Ti,Lmbd]=eig(Y(:,:,1)*Z(:,:,1)); GAM(:,:,1)=sqrtm(Lmbd); TI(:,:,1)=Ti; 
TIEIG(:,:,1)=Ti; Zcm(:,:,1)=(Ti\Y(:,:,1)/Ti.')\GAM(:,:,1);
Ups(:,:,1)=imag(s(1)).*eye(n.c)/imag(GAM(:,:,1)); 

% Inicia Processo Iterativo de Newton-Raphson
for k=2:n.Nr    
    [Ti,Lmbd,NrIte]=Mdal_EigNR(Z(:,:,k),Y(:,:,k),Ti,Lmbd,s(k),EpsTol,IterMax); 

    % Atualiaza variáveis
    TI(:,:,k)=Ti; GAM(:,:,k)=sqrtm(Lmbd); NRIte(k)=NrIte;

    % Impedância Característica Modal 
    Zcm(:,:,k) = (Ti\Y(:,:,k)/Ti.')\GAM(:,:,k); 
    
    % Velocidade de Fase Modal    
    Ups(:,:,k) = imag(s(k)).*eye(n.c)/imag(GAM(:,:,k)); 

    % Função EIG para comparação gráfica
    [TiEIG,~]=eig(Y(:,:,k)*Z(:,:,k)); TIEIG(:,:,k)=TiEIG;
end
end
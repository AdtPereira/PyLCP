function [Output,paras] = MultiConMatrix_C(x,y,a,paras,option,type,group)
%   Impedance Matrix Calculation for Round Conductors with Boundaries
%
%   Input variables
%   x           X coordinates of centre points m
%   y           Y coordinates of centre points m
%   a           Radius m
%   para.bx     Coordinates of vertical boundary m
%   para.by     Coordinates of horizental boundary m
%   para.er     Permittivity
%   para.b      Outer radius m
%   para.er2    Permittivity of insulation layer
%   para.V      Voltage V
%   para.I      Current A
%   para.Qs     Charge C
%   para.Hext	External magnetic field A/m
%   option.r0	Reference point m
%   option.Nord     Truncated order, default N = 3
%   option.k	Reflection coefficient, default k = 1
%   option.Nref Reflection time, defalt Nref = 2
%   option.a0   Boundary potential
%   type.kind   Field selection
%   type.result Type of calculation
%   type.correct Flag for boundary compensation
%   type.ExtH   Flag for external magnetic field
%
%   Output variable
%   Output:     Output based on input
%   paras       all parameters in the calculation
%
%   Reference:
%   ...
%
%   Authors:
%   Tianmingluo
%   Delft University of Technology
%
%   Email:
%   T.Luo-1@tudelft.nl
%
%   Version 1.0.238

arguments
    x  (:,1) double  % x coordinates of centre points
    y  (:,1) double  % y coordinates of centre points
    a  (:,1) double  % radius of conductors
    
    paras.bx        (1,:) double =[];
    paras.by        (1,:) double =[];
    paras.er        (1,1) double {mustBeNumeric} = 1; % permittivity of background
    paras.b         (:,1) double = [];
    paras.er2       (1,1) double {mustBeNumeric} = 1;% permittivity of layer
    paras.mur       (:,1) double = 1; % permeability of conductor
    paras.V         (:,1) double =[];
    paras.I         (:,1) double =[];
    paras.Qs        (:,1) double =[];
    paras.Hext      (1,2) double = [0,0];
    
    option.r0       double {mustBeNumeric} = 1;
    option.Nord     double {mustBeNumeric} = 3;
    option.k        double {mustBeNumeric} = 1;
    option.Nref     double {mustBeNumeric} = 2;
    option.a0       double {mustBeNumeric} = 0;
    
    type.kind       string {mustBeMember(type.kind,["C1","C2"])} = "C1";
    type.result     string {mustBeMember(type.result,["case","matrix"])} = "case";
    type.correct    logical = false;
    type.ExtH       logical = false;
    type.group      logical = false;
    
    group.gplist    (:,2) single = repmat((1:length(x))',1,2);
end
%% Constant define
e0 = 8.85e-12; % unit F/m

%% Input check
flag = 0;

warnid = 'MATLAB:nearlySingularMatrix';
warning('off',warnid);

num = length(x);

if ~isequal(length(y),num)
    flag = 1;
    errorinfo = "Coordinates must have the same length as x";
end
if ~isequal(length(a),num)
    flag = 1;
    errorinfo = "Radius must have the same length as x";
end

if type.kind == "C1"||type.kind == "C2"
    
    if type.kind == "C2"
        if isempty(paras.b)
            flag = 1;
            errorinfo = "Found empty outer radius array";
        end
    end
    if type.result == "case"
        if isempty(paras.Qs)&&isempty(paras.V)
            flag = 1;
            errorinfo = "Found empty charge or voltage array";
        end
        if ~isempty(paras.Qs)&&~isempty(paras.V)
            flag = 1;
            errorinfo = "Found double input array";
        end
        if ~isequal(length(paras.Qs),num)&&~isequal(length(paras.V),num)
            flag = 1;
            errorinfo = "Charge or voltage must have the same length as x";
        end
    end
end

if type.correct
    if isempty(paras.bx) && isempty(paras.by)
        type.correct = false;
    end
    if option.k == -1 && numel(paras.bx)<=1 && numel(paras.by)<=1
        type.correct = false;
    end
end

if type.group
    group.gplist = array2table(group.gplist);
    group.gplist.Properties.VariableNames={'NoC','NoG'};
    if  max(group.gplist.NoC)>num
        flag = 1;
        errorinfo = "the No. of conductor exceed the existing conductor";
    end
    group.gplist.NoG= categorical(group.gplist.NoG);
    gpkind = categories(group.gplist.NoG);
    gpcount = countcats(group.gplist.NoG);
    idcount = gpcount ~= 1;
    if ~nnz(idcount)
        type.group = 0;
    end
    iddel = ismember(group.gplist.NoG,gpkind(~idcount));
    group.gplist(iddel,:) = [];
    group.gplist.NoG = removecats(group.gplist.NoG,gpkind(~idcount));
    gpkind = gpkind(idcount);
end

%%
if flag == 0
    
    lnr0 = log(option.r0);
    N = option.Nord;
    k = option.k;
    Nref = option.Nref;
    a0 = option.a0;
    bx = paras.bx;
    by = paras.by;
    units = 2*N+1;
    n1 = 1:N;
    totlen = units*num;
    
    [matA0,matB0,countD] = mat_basic(x,y,bx,by,N,Nref,k);
    
    matR0 = zeros(totlen,num);
    matR0(1:units:totlen,:)= countD;
    
    if type.group % disable correction when group is used
        type.correct = false;
    end
    
    if type.correct
        [mark,bound,ind] = point_comp(x,y,bx,by,type.gap,gap.xg,gap.yg);
    else
        matB0 = matB0-matR0*lnr0;
    end
    
       
    %%
    if type.kind == "C1" || type.kind == "C2"
        b = paras.b;
        
        if type.kind == "C1"
            matf = matf_C1(a,N);
        end
        
        if type.kind == "C2"
            matf = matf_C2(a,paras.b,paras.er,paras.er2,N);
        end
        
        munit = eye(totlen);
        matA0 = matA0.*matf+munit;
        
        if type.correct
            [Acomp,Bcomp] = matComp(mark,ind,bound,x,y,matf(1,(ind-1)*units+1:ind*units,:),totlen);
            fprintf('Compensation for Boundary Potential Applied. \n');
        end
        
        idC = false(totlen,1);
        idC(1:units:totlen)=true;
        
        if type.result == "case"
            if ~isempty(paras.Qs)
                vorq = false;
                paras.V = zeros(num,1);
                if type.group
                    [matVc,matVd,matVr] = matV_C(type.kind,paras.er,paras.er2,a,b,N);
                end
            elseif ~isempty(paras.V)
                vorq = true;
                paras.Qs = zeros(num,1);
                if type.correct
                    type.correct = false;
                    fprintf('Potential correction only work with charges. \n');
                end
                [matVc,matVd,matVr] = matV_C(type.kind,paras.er,paras.er2,a,b,N);
            end
            nloop = 1;
        else
            vorq = false;
            paras.FactorAB = zeros(totlen,num);
            paras.V = zeros(num);
            paras.Qs = eye(num);
            nloop = num;
        end
        
        for idx = 1:nloop
            if ~vorq
                Qs = paras.Qs(:,idx);
                D = -Qs/2/pi/e0/paras.er;
                if type.group
                    [matA,matB] = matGroup(matA0,matB0,D,matVc,matVd-matVr*lnr0,group.gplist);
                else
                    matB = matB0*D;
                    if type.correct
                        matB = [matB;a0+D(ind)*Bcomp];
                        matR = matR0*D;
                        matA = [matA0,matR;Acomp,-D(ind)];
                    else
                        matA = matA0;
                    end
                end
            else
                matA = [matA0,-matB0;matVc,matVd-matVr*lnr0];
                matB = [zeros(totlen,1);paras.V];
            end
            
            temp_result = matA\matB;
            if vorq
                D = temp_result(totlen+1:end);
            else
                if type.correct
                    lnr0 = real(temp_result(end));
                    if ~isreal(temp_result(end))
                        temp_result(1:totlen) = temp_result(1:totlen)+1i*imag(temp_result(end))*matR;
                    end
                end
                if type.group
                    D(group.gplist.NoC) = temp_result(totlen+1:end);
                end
            end
            paras.FactorAB(:,idx) = temp_result(1:totlen);
            option.r0(1,idx) = exp(lnr0);
            C = paras.FactorAB(idC,idx);
            
            if vorq || type.group
                    paras.Qs(:,idx) = -D*2*pi*e0*paras.er;
            end
            if ~vorq
                if type.kind == "C1"
                    paras.V(:,idx) = C+D.*(log(a)-lnr0);
                else
                    paras.V(:,idx) = C+D.*(log(b)-lnr0+paras.er/paras.er2*log(a./b));
                end
            end
        end
        paras.I = 1j*paras.Qs*w;
        
        S = paras.V.*conj(paras.I);
        paras.P = real(S);
        paras.Q = imag(S);
        if type.result == "matrix"
            Cap = inv(paras.V);
        end
    end
    %%

    if type.result == "matrix"
        if type.kind == "L"
            Output = Imp;
        end
        if type.kind == "C1" || type.kind == "C2"
            Output = Cap;
        end
    else
        if type.kind == "L"
            if exist('Imp','var')
                Output = Imp;
            else
                Output = paras.P;
            end
        else
            if vorq
                Output = paras.Qs;
            else
                Output = paras.V;
            end
        end
    end
    
    paras.x = x;
    paras.y = y;
    paras.a = a;
    paras.option = option;
    paras.type = type;
    paras.gap = gap;
else
    error(errorinfo);
end

end

function matf = matf_C1(a,N)
%--- Matrix converting A'' and B'' coefficients to A' and B' in C1 situation ---

num = length(a);
units = 2*N+1;

matf = cell(num,num);
[matn,~] = meshgrid((1:N)',(0:N)');
um = ones(units,1);

for id = 1:num
    pf = -a(id).^(2*matn);
    pf2 = pf(2:end,:);
    matf(:,id)= {[um,[repmat(pf,1,2);repmat(pf2,1,2)]]};
end

matf = cell2mat(matf);
end

function matf = matf_C2(a,b,er,er2,N)
%--- Matrix converting A'' and B'' coefficients to A' and B' in C2 situation ---

num = length(a);
units = 2*N+1;

matf = cell(num,num);
[matn,~] = meshgrid((1:N)',(0:N)');
um = ones(units,1);

erp = er+er2;
erm = er2-er;

for id = 1:num
    part1 = a(id).^(2*matn);
    part2 = b(id).^(2*matn);
    pf = -part2.*(erm*part2+erp*part1)./(erp*part2+erm*part1);
    pf2 = pf(2:end,:);
    matf(:,id)= {[um,[repmat(pf,1,2);repmat(pf2,1,2)]]};
end

matf = cell2mat(matf);
end

function [matVc,matVd,matVr] = matV_C(type,er,er2,a,b,N)
arguments
    type string {mustBeMember(type,["C1","C2"])} = "C1";
    er = 1;
    er2 = 1;
    a double = [];
    b double = [];
    N double = 3;
end
% e0 = 8.85e-12; % unit F/m
num = length(a);
units = 2*N+1;

matVc = cell(num,num);
matVc(:) = {zeros(1,units)};
for k = 1:num
    matVc(k,k) = {eye(1,units)};
end
matVc = cell2mat(matVc);

switch type
    case "C1"
        chi = log(a);
    case "C2"
        chi = log(b)+er/er2*log(a./b);
end

matVd = diag(chi);
matVr = eye(num);
end

function [mark,bound,ind] = point_comp(x,y,bx,by,flag,xg,yg)
%--- Select cell for compensation ---
arguments
    x
    y
    bx      (1,:) = [];
    by      (1,:) = [];
    flag    = false;
    xg      (1,:) = [];
    yg      (1,:) = [];
end
if ~isempty(bx)
    dx = abs(x-bx);
    valx = min(dx,[],'all');
    midx = mean(bx);
end
if ~isempty(by)
    dy = abs(y-by);
    valy = min(dy,[],'all');
    midy = mean(by);
end
if flag
    dxgap = min(abs(x-xg),[],2);
    dygap = min(abs(y-yg),[],2);
end

if ~isempty(bx)&&~isempty(by)
    if valx<valy
        mark = 'x';
        if ~flag
            [~,ind,bound_id] = process(dx,valx,abs(y-midy),1);
        else
            [~,ind,bound_id] = process(dx,valx,dygap,0);
        end
        bound = bx(bound_id);
    elseif valx>valy
        mark = 'y';
        if ~flag
            [~,ind,bound_id] = process(dy,valy,abs(x-midx),1);
        else
            [~,ind,bound_id] = process(dy,valy,dxgap,0);
        end
        bound = by(bound_id);
    else
        note = 0;
        if ~flag
            [valx2,indx,bound_idx] = process(dx,valx,abs(y-midy),1);
            [valy2,indy,bound_idy] = process(dy,valy,abs(x-midx),1);
            if valx2 <= valy2
                note = 1;
            end
        else
            [valx2,indx,bound_idx] = process(dx,valx,dygap,0);
            [valy2,indy,bound_idy] = process(dy,valy,dxgap,0);
            if valx2 >= valy2
                note = 1;
            end
        end
        if note
            mark = 'x';
            ind = indx;
            bound = bx(bound_idx);
        else
            mark = 'y';
            ind = indy;
            bound = by(bound_idy);
        end
    end
elseif ~isempty(bx)
    mark = 'x';
    if ~flag
        [~,ind,bound_id] = process(dx,valx,abs(y-mean(y)),1);
    else
        [~,ind,bound_id] = process(dx,valx,dygap,0);
    end
    bound = bx(bound_id);
elseif ~isempty(by)
    mark = 'y';
    if ~flag
        [~,ind,bound_id] = process(dy,valy,abs(x-mean(x)),1);
    else
        [~,ind,bound_id] = process(dy,valy,dxgap,0);
    end
    bound = by(bound_id);
end
    function [val,C_id,bound_id] = process(a,valin,clue,type)
        ind_1 = find(a==valin);
        [row,col] = ind2sub(size(a),ind_1);
        if type
            [val,ind_2] = min(clue(row));
        else
            [val,ind_2] = max(clue(row));
        end
        C_id = row(ind_2);
        bound_id = col(ind_2);
    end
end

function [matcomp,Bcomp] = matComp(mark,ind,bound,x,y,matf,totlen)
%--- parameters in compensation equation ---
[~,units,num2] = size(matf);
N =(units-1)/2;
n1 = 1:N;
matcomp = zeros(1,totlen,num2);

switch mark
    case 'x'
        xm = bound-x(ind);
        A21p1 = xm.^(n1);
        A21p2 = xm.^(-n1);
        Bcomp = -log(abs(xm));
        
    case 'y'
        ym = bound-y(ind);
        A21p1 = (1i*ym).^(n1);
        A21p2 = (-1i*ym).^(-n1);
        Bcomp = -log(abs(ym));
end

A21p1 = [0,real(A21p1),imag(A21p1)];
A21p2 = [1,real(A21p2),imag(A21p2)];

matcomp(1,(ind-1)*units+1:ind*units,:) = A21p1+A21p2.*matf;
end

function [matA, matB] = matGroup(matA0,matB0,D,matVc,matVd,gplist)
num = length(D);
[totlen,~,num2] = size(matA0);

numgp = numel(gplist.NoC);

gpkind = categories(gplist.NoG);
nog = numel(gpkind);

matA = [matA0,-repmat(matB0(:,gplist.NoC),1,1,num2)];
matA = [matA;zeros(numgp,totlen+numgp,num2)];

matB0(:,gplist.NoC) = [];

tempD = D;
tempD(gplist.NoC) = [];

matB = matB0*tempD;
matB = [matB;zeros(numgp,1)];

id = totlen+1;

for k = 1:nog
    kind = gpkind(k);
    list = ismember(gplist.NoG,kind);
    matB(id) = sum(D(list));
    rowD = zeros(1,numgp);
    rowD(list) = 1;
    matA(id,totlen+1:end,:) = repmat(rowD,1,1,num2);
    
    idx = find(list == 1);
    rowVc = matVc(idx(2:end),:,:) - matVc(idx(1),:,:);
    rowVd = matVd(idx(2:end),gplist.NoC,:) - matVd(idx(1),gplist.NoC,:);
    
    matA(id+1:id+nnz(list)-1,:,:) = [rowVc,rowVd];
    id = id + nnz(list);
end

end
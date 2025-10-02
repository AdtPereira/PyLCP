function [matA,matB,countD] =  mat_basic(x,y,bx,by,N,Nref,k)
%   Matrix Geometric Parts Construction for 2D Round Conductors Approach
%
%   Input variables
%   N           Truncated order, default N = 3
%   x           X coordinates of centre points m
%   y           Y coordinates of centre points m
%   bx          Coordinates of vertical boundary m
%   by          Coordinates of horizental boundary m
%   Nref        Number of reflection
%   k           Reflection coefficient, default k = 1
%
%   Output variable
%   matA:       Geometric part on the left side
%   matB:       Geometric part on the right side
%   countD:     Contribution times of each round cell
%
%   Reference:
%   ...
%
%   Authors:
%   Tianmingluo
%   Delft University of Technology
%
%   Email:
%   T.Luo-1@tudelft.n

arguments
    x       (:,1) double  % x coordinates of centre points
    y       (:,1) double  % y coordinates of centre points
    bx      (1,:) double
    by      (1,:) double
    N       (1,1) double {mustBeNumeric};
    Nref    (1,1) double {mustBeNumeric}
    k       (1,1) double {mustBeNumeric}
end

if nnz(isnan(bx))
    id = isnan(bx);
    bx(id) = [];
end

if nnz(isnan(by))
    id = isnan(by);
    by(id) = [];
end

units = 2*N+1;

num = length(x);

matAij = cell(num,num);

matAij(:) = {zeros(units)};

matB = cell(num,num);

matB(:) = {zeros(units,1)};

countD = zeros(num);

[n,m] = meshgrid((1:N)',(0:N)');

n1 = (1:N)';

factor = factorial(n+m-1)./factorial(n-1)./factorial(m);

antiB = [1;(-1).^n1;(-1).^n1];

nm = -(n+m);

antiA = (-1).^nm;

antiA2 = antiA(2:end,:);

geop1 = (-1).^n;

umij = zeros(2*N+1,1);

nx = length(bx);
ny = length(by);
%%
for idx = 1:num
    
    dx = x-x(idx);
    dy = y-y(idx);
    
    iterate_zero(dx,dy,idx);
    
end

matA = matAij;
%%
if (nx+ny)~=0
    
    if nx == 2
        w = bx(2)-bx(1);
        w(1,1,2) = -w(1,1,1);
    end
    
    if ny == 2
        h = by(2)-by(1);
        h(1,1,2) = -h(1,1,1);
    end
    
    matAij = cell(num,num);
    matAij(:) = {zeros(units)};
    
    for idx = 1:num
        if nx~= 0
            dx1 = 2*bx(1)-x-x(idx);
            if nx == 2
                dx2 = 2*bx(2)-x-x(idx);
            end
        end
        if ny~=0
            dy1 = 2*by(1)-y-y(idx);
            if ny ==2
                dy2 = 2*by(2)-y-y(idx);
            end
        end
        dy = y-y(idx);
        dx = x-x(idx);
        
        if nx~=0
            iterate_one(dx1,dy,"x",idx,k^1);
            if nx == 2
                iterate_one(dx2,dy,"x",idx,k^1);
                iterate_two(dx,dx2,dy,w,"xnone",idx,k);
            end
        end
        
        if ny~=0
            iterate_one(dx,dy1,"y",idx,k^1);
            if ny == 2
                iterate_one(dx,dy2,"y",idx,k^1);
                iterate_two(dy,dy2,dx,h,"ynone",idx,k);
            end
        end
        
        if ny==1 && nx ==1
            iterate_one(dx1,dy1,"xy",idx,k^2);
        end
        
        if ny ==1 && nx ==2
            iterate_two(dx,dx2,dy1,w,"yxy",idx,k);
        end
        
        if nx == 1 && ny ==2
            iterate_two(dy,dy2,dx1,h,"xxy",idx,k);
        end
        
        if nx+ny == 4
            iterate_four(dx,dx2,dy,dy2,w,h,idx,k);
        end
    end
    
    matAxy = matAij;
end
%%
matA = cell2mat(matA);

if (nx+ny)~=0
    matA = matA+ cell2mat(matAxy);
end

matB = cell2mat(matB);
%%
    function iterate_zero(dx,dy,idx)
        for idx2 = idx:num
            
            xin = dx(idx2);
            yin = dy(idx2);
            
            [Bij,Bji,Aij,Aji,Dij] = computblock(xin,yin,1,"none");
            
            if idx ~= idx2
                adding(idx,idx2,Aij,Bij,Dij);
                adding(idx2,idx,Aji,Bji,Dij);
            end
        end
    end

    function iterate_one(dx,dy,type,idx,k0)
        for idx2 = idx:num
            
            xin = dx(idx2);
            yin = dy(idx2);
            
            [Bij,Bji,Aij,Aji,Dij] = computblock(xin,yin,k0,type);
            
            adding(idx,idx2,Aij,Bij,Dij);
            if idx ~= idx2
                adding(idx2,idx,Aji,Bji,Dij);
            end
        end
    end

    function iterate_two(dx,dx2,dy,w,type,idx,k0)
        for idx2 = idx:num
            yin = dy(idx2);
            for time = 2:Nref
                if ~mod(time,2)
                    if contains(type,"none")
                        xin = dx(idx2)+time*w;
                    else
                        vmove = reshape([time-2,time],[1,1,2]);
                        xin = dx2(idx2)+vmove.*w;
                    end
                    
                    switch type
                        case "xnone"
                            [Bij,Bji,Aij,Aji,Dij] = computblock(xin,yin,k0^time,"none");
                        case "ynone"
                            [Bij,Bji,Aij,Aji,Dij] = computblock(yin,xin,k0^time,"none");
                        case "yxy"
                            [Bij,Bji,Aij,Aji,Dij] = computblock(xin,yin,k0^time,"xy");
                        case "xxy"
                            [Bij,Bji,Aij,Aji,Dij] = computblock(yin,xin,k0^time,"xy");
                    end
                else
                    if contains(type,"none")
                        vmove = reshape([time-1,time+1],[1,1,2]);
                        xin = dx2(idx2)+vmove.*w;
                    else
                        xin = dx(idx2)+(time-1)*w;
                    end
                    
                    switch type
                        case "xnone"
                            [Bij,Bji,Aij,Aji,Dij] = computblock(xin,yin,k0^time,"x");
                        case "ynone"
                            [Bij,Bji,Aij,Aji,Dij] = computblock(yin,xin,k0^time,"y");
                        case "yxy"
                            [Bij,Bji,Aij,Aji,Dij] = computblock(xin,yin,k0^time,"y");
                        case "xxy"
                            [Bij,Bji,Aij,Aji,Dij] = computblock(yin,xin,k0^time,"x");
                    end
                end
                
                adding(idx,idx2,Aij,Bij,Dij);
                if idx ~= idx2
                    adding(idx2,idx,Aji,Bji,Dij);
                end
            end
        end
    end

    function iterate_four(dx,dx2,dy,dy2,w,h,idx,k0)
        for idx2 = idx:num
            for time = 2:Nref
                
                mx = 1:time-1;
                my = time-mx;
                
                for idx3 = 1:time-1
                    
                    if mod(mx(idx3),2) == 1 && mod(my(idx3),2) == 1
                        
                        xmove = reshape([mx(idx3)-1,mx(idx3)+1],[1,1,2]);
                        ymove = reshape([my(idx3)-1,my(idx3)+1],[1,1,2]);
                        
                        xin = dx2(idx2)+repmat(xmove,1,1,2).*repmat(w,1,1,2);
                        yin = dy2(idx2)+reshape(repmat(ymove,1,2).*repmat(h,1,2),1,1,4);
                        [Bij,Bji,Aij,Aji,Dij] = computblock(xin,yin,k0^time,"xy");
                        
                    elseif mod(mx(idx3),2) == 0 && mod(my(idx3),2) == 1
                        
                        ymove = reshape([my(idx3)-1,my(idx3)+1],[1,1,2]);
                        
                        xin = dx(idx2)+ mx(idx3)*repmat(w,1,1,2);
                        yin = dy2(idx2)+reshape(repmat(ymove,1,2).*repmat(h,1,2),1,1,4);
                        [Bij,Bji,Aij,Aji,Dij] = computblock(xin,yin,k0^time,"y");
                        
                    elseif mod(mx(idx3),2) == 1 && mod(my(idx3),2) == 0
                        
                        xmove = reshape([mx(idx3)-1,mx(idx3)+1],[1,1,2]);
                        
                        xin = dx2(idx2)+repmat(xmove,1,1,2).*repmat(w,1,1,2);
                        yin = dy(idx2)+reshape(my(idx3)*repmat(h,1,2),1,1,4);
                        [Bij,Bji,Aij,Aji,Dij] = computblock(xin,yin,k0^time,"x");
                        
                    else
                        xin = dx(idx2)+ mx(idx3)*repmat(w,1,1,2);
                        yin = dy(idx2)+reshape(my(idx3)*repmat(h,1,2),1,1,4);
                        [Bij,Bji,Aij,Aji,Dij] = computblock(xin,yin,k0^time,"none");
                    end
                    
                    adding(idx,idx2,Aij,Bij,Dij);
                    if idx ~= idx2
                        adding(idx2,idx,Aji,Bji,Dij);
                    end
                end
            end
        end
    end

    function iterate_four_v2(dx,dx2,dy,dy2,w,h,idx,k0)
        for idx2 = idx:num
            for idx3 = 1:Nref
                for idx4 = 1:Nref
                    if mod(idx3,2) == 1 && mod(idx4,2) == 1
                        
                        xmove = reshape([idx3-1,idx3+1],[1,1,2]);
                        ymove = reshape([idx4-1,idx4+1],[1,1,2]);
                        
                        xin = dx2(idx2)+repmat(xmove,1,1,2).*repmat(w,1,1,2);
                        yin = dy2(idx2)+reshape(repmat(ymove,1,2).*repmat(h,1,2),1,1,4);
                        [Bij,Bji,Aij,Aji,Dij] = computblock(xin,yin,k0^(idx3+idx4),"xy");
                        
                    elseif mod(idx3,2) == 0 && mod(idx4,2) == 1
                        
                        ymove = reshape([idx4-1,idx4+1],[1,1,2]);
                        
                        xin = dx(idx2)+ idx3*repmat(w,1,1,2);
                        yin = dy2(idx2)+reshape(repmat(ymove,1,2).*repmat(h,1,2),1,1,4);
                        [Bij,Bji,Aij,Aji,Dij] = computblock(xin,yin,k0^(idx3+idx4),"y");
                        
                    elseif mod(idx3,2) == 1 && mod(idx4,2) == 0
                        
                        xmove = reshape([idx3-1,idx3+1],[1,1,2]);
                        
                        xin = dx2(idx2)+repmat(xmove,1,1,2).*repmat(w,1,1,2);
                        yin = dy(idx2)+reshape(idx4*repmat(h,1,2),1,1,4);
                        [Bij,Bji,Aij,Aji,Dij] = computblock(xin,yin,k0^(idx3+idx4),"x");
                        
                    else
                        xin = dx(idx2)+ idx3*repmat(w,1,1,2);
                        yin = dy(idx2)+reshape(idx4*repmat(h,1,2),1,1,4);
                        [Bij,Bji,Aij,Aji,Dij] = computblock(xin,yin,k0^(idx3+idx4),"none");
                    end
                    
                    adding(idx,idx2,Aij,Bij,Dij);
                    if idx ~= idx2
                        adding(idx2,idx,Aji,Bji,Dij);
                    end
                end
            end
        end
    end

    function adding(i,j,A,B,D)
        matAij{i,j} = matAij{i,j}+A;
        matB{i,j} = matB{i,j}+B;
        countD(i,j) = countD(i,j)+D;
    end

    function [Bij,Bji,Aij,Aji,Dij] = computblock(xin,yin,k,type)
        
        ninput = max(length(xin),length(yin));
        
        [B,rfgeo,ifgeo] = Fgeo(xin,yin,n1);
        
        Dij = k*ninput;
        
        Bij = -k*[B;rfgeo;ifgeo];
        
        [Aij,Aji] = contri_ij(xin,yin,k,type);
        
        switch type
            case "none"
                Bji = Bij.*antiB;
            case "x"
                Bji = -k*[B;rfgeo;-ifgeo];
            case "y"
                Bji = -k*[B;rfgeo;-ifgeo].*antiB;
            case "xy"
                Bji = Bij;
        end
    end

    function [Mij,Mji] = contri_ij(dx,dy,k,type)
        
        % contribution from Aj to Ai, N is order number
        % total matrix size is (2N+1)*(2N+1)
        % Matrix looks like | Cj2Ci   Ajn2C   Bjn2C   |
        %                   | Cj2Ain  Ajn2Ain Bjn2Ain |
        %                   | Cj2Bin  Ajn2Bin Bjn2Bin |
        switch type
            case "none"
                pgeo = sum(k.*geop1.*(dx-1i*dy).^nm.*factor,3);
            case "x"
                pgeo = sum(k.*(dx+1i*dy).^nm.*factor,3);
            case "y"
                pgeo = sum(k.*geop1.*(dx+1i*dy).^nm.*factor,3);
            case "xy"
                pgeo = sum(k.*(dx-1i*dy).^nm.*factor,3);
        end
        
        pgeo2 = pgeo(2:end,:);
        
        mr = real(pgeo);
        mrsub = real(pgeo2);
        
        mi = imag(pgeo);
        misub = imag(pgeo2);
        
        if type == "none" || type == "y"
            mr2 = mr.*antiA;
            mrsub2 = mrsub.*antiA2;
            mi2 = mi.*antiA;
            misub2 = misub.*antiA2;
        else
            mr2 = mr;
            mrsub2 = mrsub;
            mi2 = mi;
            misub2 = misub;
        end
        
        switch type
            case "none"
                Mij = [umij,[-mr,-mi;-misub,mrsub]];
                Mji = [umij,[-mr2,-mi2;-misub2,mrsub2]];
            case "x"
                Mij = [umij,[-mr,-mi; misub,-mrsub]];
                Mji = [umij,[-mr2,mi2; -misub2, -mrsub2]];
            case "y"
                Mij = [umij,[-mr,-mi; misub, -mrsub]];
                Mji = [umij,[-mr2,mi2; -misub2, -mrsub2]];
            case "xy"
                Mij = [umij,[-mr,-mi;-misub,mrsub]];
                Mji = [umij,[-mr2,-mi2;-misub2,mrsub2]];
        end
    end
end

function [B0,fgeor,fgeoi] = Fgeo(x,y,n)

arguments
    x  (1,:) double
    y  (1,:) double
    n  (:,1) double
end

Fac = sum((x-1i*y).^-n./n,2);

fgeor = real(Fac);
fgeoi = imag(Fac);

B0 = sum(-1/2*log(x.^2+y.^2));

end

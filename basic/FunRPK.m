function [Fs] = FunRPK(r,p,k,s)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

for a=1:length(s)
    for n=1:length(r)
        Zn(n) = r(n)/(s(a)-p(n));
    end
    Fs(a) = k + sum(Zn);
end
end
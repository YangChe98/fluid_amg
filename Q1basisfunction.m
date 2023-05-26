function [phi,phix,phiy]=Q1basisfunction(n,xi,eta)
%% n basisfunction in Q1
if n==1
    phi=(1-xi).*(1-eta)/4;
    phix=-(1-eta)/4;
    phiy=-(1-xi)/4;
elseif n==2
    phi=(1+xi).*(1-eta)/4;
    phix=(1-eta)/4;
    phiy=-(1+xi)/4;
elseif n==3
    phi=(1+xi).*(1+eta)/4;
    phix=(1+eta)/4;
    phiy=(1+xi)/4;
else
    phi=(1-xi).*(1+eta)/4;
    phix=-(1+eta)/4;
    phiy=(1-xi)/4;
end
end
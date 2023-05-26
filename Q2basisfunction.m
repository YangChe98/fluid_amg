function [phi,phix,phiy]=Q2basisfunction(n,xi,eta)


%% n basisfunction in Q2
if n==1
    phi=(1-xi).*(1-eta).*eta.*xi/4;
    phix=(eta.*xi.*(eta - 1))/4 + (eta.*(eta - 1).*(xi - 1))/4;
    phiy=(eta.*xi.*(xi - 1))/4 + (xi.*(eta - 1).*(xi - 1))/4;
elseif n==2
    phi=-(1+xi).*(1-eta).*eta.*xi/4;
    phix=(eta.*xi.*(eta - 1))/4 + (eta.*(eta - 1).*(xi + 1))/4;
    phiy=(eta.*xi.*(xi + 1))/4 + (xi.*(eta - 1).*(xi + 1))/4;
elseif n==3
    phi=(1+xi).*(1+eta).*eta.*xi/4;
    phix=(eta.*xi.*(eta + 1))/4 + (eta.*(eta + 1).*(xi + 1))/4;
    phiy=(eta.*xi.*(xi + 1))/4 + (xi.*(eta + 1).*(xi + 1))/4;
elseif n==4
    phi=-(1-xi).*(1+eta).*eta.*xi/4;
    phix=(eta.*xi.*(eta + 1))/4 + (eta.*(eta + 1).*(xi - 1))/4;
    phiy=(eta.*xi.*(xi - 1))/4 + (xi.*(eta + 1).*(xi - 1))/4;
 elseif n==5
    phi=-(1-xi.^2).*(1-eta).*eta/2;
    phix=-eta.*xi.*(eta - 1);
    phiy=-((xi.^2 - 1).*(eta - 1))/2 - (eta.*(xi.^2 - 1))/2;
  elseif n==6
    phi=(1-eta.^2).*(1+xi).*xi/2;
    phix=- ((eta.^2 - 1).*(xi + 1))/2 - (xi.*(eta.^2 - 1))/2;
    phiy=-eta.*xi.*(xi + 1);
  elseif n==7
    phi=(1-xi.^2).*(1+eta).*eta/2;
    phix=-eta.*xi.*(eta + 1);
    phiy=- ((xi.^2 - 1).*(eta + 1))/2 - (eta.*(xi.^2 - 1))/2;
  elseif n==8
    phi=-(1-eta.^2).*(1-xi).*xi/2;
    phix=- ((eta.^2 - 1).*(xi - 1))/2 - (xi.*(eta.^2 - 1))/2;
    phiy=-eta.*xi.*(xi - 1);
else
     phi=(1-eta.^2).*(1-xi.^2);
    phix=2.*xi.*(eta.^2 - 1);
    phiy=2.*eta.*(xi.^2 - 1);
end
end
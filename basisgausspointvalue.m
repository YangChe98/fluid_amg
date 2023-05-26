function [phigaussvalue]=basisgausspointvalue()

syms x y;

gausspoint=[-sqrt(3/5),0,sqrt(3/5)];

xgausspoint2d=repmat(gausspoint,3,1);
ygausspoint2d=repmat(gausspoint.',1,3);
xgausspoint2d=reshape(xgausspoint2d,1,[]);
ygausspoint2d=reshape(ygausspoint2d,1,[]);


  phigaussvalue=zeros(9,9);
for i=1:9
     phigaussvalue(i,:)=Q2basisfunction(i,xgausspoint2d,ygausspoint2d);


end

save phigaussvalue.mat phigaussvalue;

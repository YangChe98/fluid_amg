function plotuvp0(x_min,x_n,x_max,y_min,y_n,y_max,basisfunctionweight,node_number,ubasis_function_number,pbasis_function_number)



xplot=x_min:(1/(2*x_n)):x_max;
 yplot=y_min:(1/(2*y_n)):y_max;
xplot1=x_min:(1/(x_n)):(x_max);
 yplot1=y_min:(1/(y_n)):(y_max);


[xplot,yplot]=meshgrid(xplot,yplot);
[xplot1,yplot1]=meshgrid(xplot1,yplot1);
% xplot1=xplot1.';
% yplot1=yplot1.';
% xplot=xplot.';
% yplot=yplot.';

uvalue=basisfunctionweight(1:node_number,1);

uvalue=reshape(uvalue,size(xplot.'));
vvalue=basisfunctionweight(node_number+1:node_number*2,1);
vvalue=reshape(vvalue,size(xplot.'));

pvalue=basisfunctionweight(node_number*2+1:ubasis_function_number+pbasis_function_number,1);
 pvalue=reshape(pvalue,size(xplot1.'));
figure(1)

streamslice(xplot,yplot,uvalue.',vvalue.')

figure(2)

mesh(xplot1,yplot1,pvalue.')
colorbar
figure(3)

mesh(xplot,yplot,uvalue.')
colorbar
figure(4)

mesh(xplot,yplot,vvalue.')
colorbar
% uvalueana=function_u(xplot,yplot);
% vvalueana=function_v(xplot,yplot);
% pvalueana=function_p(xplot1,yplot1);
% figure(5)
% 
% mesh(xplot1,yplot1,pvalueana)
% colorbar
% figure(6)
% 
% mesh(xplot,yplot,uvalueana)
% colorbar
% figure(7)
% 
% mesh(xplot,yplot,vvalueana)
% figure(8)
% streamslice(xplot,yplot,uvalueana,vvalueana)

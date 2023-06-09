function plotuvp1(x_min,x_n,x_max,y_min,y_n,y_max,basisfunctionweight,node_number,ubasis_function_number,pbasis_function_number)



xplot=x_min:(1/(2*x_n)):x_max;
 yplot=y_min:(1/(2*y_n)):y_max;
xplot1=x_min:(1/(x_n)):(x_max);
 yplot1=y_min:(1/(y_n)):(y_max);

[xplot,yplot]=meshgrid(xplot,yplot);
[xplot1,yplot1]=meshgrid(xplot1,yplot1);

uvalue=basisfunctionweight(1:node_number,1);
uvalue=reshape(uvalue,size(xplot.'));
vvalue=basisfunctionweight(node_number+1:node_number*2,1);
vvalue=reshape(vvalue,size(xplot.'));

pvalue=basisfunctionweight(node_number*2+1:ubasis_function_number+pbasis_function_number,1);
 pvalue=reshape(pvalue,size(xplot1.'));
figure(1)

streamslice(xplot,yplot,uvalue.',vvalue.')

figure(2)
%surf(xplot1,yplot1,pvalue.')
%mesh(xplot1,yplot1,pvalue.')
%shading interp
contourf(xplot,yplot,sqrt(uvalue.'.^2+vvalue.'.^2))
colorbar
figure(3)

mesh(xplot,yplot,uvalue.')
shading interp;
colorbar
figure(4)

mesh(xplot,yplot,vvalue.')
shading interp;
colorbar
uvalueana=function_u(xplot,yplot);
vvalueana=function_v(xplot,yplot);
pvalueana=function_p(xplot1,yplot1);
figure(5)

contourf(xplot,yplot,sqrt(uvalueana.^2+vvalueana.^2))
shading interp;
colorbar
figure(6)

mesh(xplot,yplot,uvalueana)
shading interp;
colorbar
figure(7)

mesh(xplot,yplot,vvalueana)
shading interp;
colorbar
figure(8)
streamslice(xplot,yplot,uvalueana,vvalueana)
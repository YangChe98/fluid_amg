
%% 2D Quadrilateral Finite Element Methods Using Discontinuous Pressures
% u velocity using Q2 pressure using Q1 ,[X,Y] is the mesh of u, [X0, Y0]
% is the mesh of p. 0 is for p. without 0 is for u.

clc
clear
tic
nu=1;
syms xi eta 

x_max=1;
x_min=0;
y_max=1;
y_min=0;
x_n=2^4;
y_n=2^4;
localbasisfunctionnumber1=9;
localbasisfunctionnumber2=4;
node_number=(2*x_n+1)*(2*y_n+1);
element_number=x_n*y_n;
ubasis_function_number=node_number*2;
pbasis_function_number=(x_n+1)*(y_n+1);
[X0,Y0,node_coordinate0,element_coordinate0,X,Y,node_coordinate,element_coordinate]=meshgenerate(x_min,x_max,y_min,y_max,x_n,y_n);
%[x,res]=solution(nu,x_n,y_n,element_number,localbasisfunctionnumber1,element_coordinate0,element_coordinate,node_coordinate,ubasis_function_number,pbasis_function_number,node_number);
[x,res]=solution6(nu,x_n,y_n,element_number,localbasisfunctionnumber1,element_coordinate0,element_coordinate,node_coordinate,ubasis_function_number,pbasis_function_number,node_number);
%[x,res]=solution2(nu,x_n,y_n,element_number,localbasisfunctionnumber1,element_coordinate0,element_coordinate,node_coordinate,ubasis_function_number,pbasis_function_number,node_number);
%[x,res]=solution3(nu,x_n,y_n,element_number,localbasisfunctionnumber1,element_coordinate0,element_coordinate,node_coordinate,ubasis_function_number,pbasis_function_number,node_number);
%[x,res]=solution5(nu,x_n,y_n,element_number,localbasisfunctionnumber1,element_coordinate0,element_coordinate,node_coordinate,ubasis_function_number,pbasis_function_number,node_number);
%Ac=RScoarsen(Atotal);
%Ac=Beckcoarsen(Atotal);
%basisfunctionweight=solution(nu,x_n,y_n,element_number,localbasisfunctionnumber1,element_coordinate0,element_coordinate,node_coordinate,ubasis_function_number,pbasis_function_number,node_number);
%plotuvp(x_min,x_n,x_max,y_min,y_n,y_max,basisfunctionweight,node_number,ubasis_function_number,pbasis_function_number)
 %plotuvp0(x_min,x_n,x_max,y_min,y_n,y_max,basisfunctionweight,node_number,ubasis_function_number,pbasis_function_number)
toc
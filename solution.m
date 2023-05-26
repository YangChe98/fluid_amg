function [x,res]=solution(nu,x_n,y_n,element_number,localbasisfunctionnumber1,element_coordinate0,element_coordinate,node_coordinate,ubasis_function_number,pbasis_function_number,node_number)


load A11reference.mat;
load A12reference.mat;
load A22reference.mat;

load B11reference.mat;
load B12reference.mat;





load phigaussvalue.mat;


%%%%%%%%%%%%%%% gauss quadrature %%%%%%%%%%
gaussweight=[5/9,8/9,5/9];
gausspoint=[-sqrt(3/5),0,sqrt(3/5)];

xgausspoint2d=repmat(gausspoint,3,1);
ygausspoint2d=repmat(gausspoint.',1,3);
xgausspoint2d=reshape(xgausspoint2d,1,[]);
ygausspoint2d=reshape(ygausspoint2d,1,[]);

gaussweight2d=gaussweight.*gaussweight.';
gaussweight2d=reshape(gaussweight2d,1,[]);
gaussweight2d=repmat(gaussweight2d,localbasisfunctionnumber1,1);

A=sparse(ubasis_function_number,ubasis_function_number);
B=sparse(ubasis_function_number,pbasis_function_number);
right=sparse(ubasis_function_number+pbasis_function_number,1);
imatrix=[];


J=[1/(x_n*2),0;0,1/(y_n*2)];
 detJ=det(J);
 detJ=abs(detJ);
 invJ=inv(J);
 invJT=invJ.';
 C=invJ*invJT;
 Alocal=nu*detJ*(C(1,1)*(A11reference)+C(1,2)*(A12reference+A12reference.')+C(2,2)*(A22reference));
    B1local=detJ*(invJ(1,1)*B11reference+invJ(2,1)*B12reference);
      B2local=detJ*(invJ(1,2)*B11reference+invJ(2,2)*B12reference);
%%%%%%%%%%%%gengerate  matrix A   B right
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i=1:element_number 
    local_coordinate1=node_coordinate(element_coordinate(i,1),1:2);
    local_coordinate2=node_coordinate(element_coordinate(i,2),1:2);
    local_coordinate3=node_coordinate(element_coordinate(i,3),1:2);
    local_coordinate4=node_coordinate(element_coordinate(i,4),1:2);
    
    Fkb=(local_coordinate1+local_coordinate3).'/2;
   
    xlocalgausspoint=J(1,1)*xgausspoint2d+J(1,2)*ygausspoint2d+Fkb(1,1)*ones(size(xgausspoint2d));
    ylocalgausspoint=J(2,1)*xgausspoint2d+J(2,2)*ygausspoint2d+Fkb(2,1)*ones(size(xgausspoint2d));
  
    fvalue=f(xlocalgausspoint,ylocalgausspoint);

    
     
   
    rightlocal1=detJ*(fvalue(1,:).*phigaussvalue).*gaussweight2d;
    rightlocal2=detJ*(fvalue(2,:).*phigaussvalue).*gaussweight2d;
    rightlocal1=sum(rightlocal1,2);
    rightlocal2=sum(rightlocal2,2);
    
 
A(element_coordinate(i,:),element_coordinate(i,:))=A(element_coordinate(i,:),element_coordinate(i,:))+Alocal;
 A(element_coordinate(i,:)+node_number,element_coordinate(i,:)+node_number)=A(element_coordinate(i,:)+node_number,element_coordinate(i,:)+node_number)+Alocal;
 B(element_coordinate(i,:),element_coordinate0(i,:))=B(element_coordinate(i,:),element_coordinate0(i,:))+B1local;
 B(element_coordinate(i,:)+node_number,element_coordinate0(i,:))=B(element_coordinate(i,:)+node_number,element_coordinate0(i,:))+B2local;
  right(element_coordinate(i,:),1) =right(element_coordinate(i,:),1)+rightlocal1;
   right(element_coordinate(i,:)+node_number,1) =right(element_coordinate(i,:)+node_number,1)+rightlocal2;
end


Atotal=[A,B;B.',zeros(pbasis_function_number,pbasis_function_number)];

%%%%%%%%%%%%%%%%%%%%% boundary condition %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i=1:node_number
    if(node_coordinate(i,3)==1)
        right=right-Atotal(:,i)*node_coordinate(i,4)-Atotal(:,node_number+i)*node_coordinate(i,5);
        Atotal(:,i)=0;
        Atotal(i,:)=0;
        
        Atotal(:,node_number+i)=0;
        Atotal(node_number+i,:)=0;
        Atotal(i,i)=1;
        Atotal(node_number+i,node_number+i)=1;
        right(i,1)=node_coordinate(i,4);
        right(node_number+i,1)=node_coordinate(i,5);
    end

end

 right=right-Atotal(:,ubasis_function_number+1)*node_coordinate(1,6);
        Atotal(:,ubasis_function_number+1)=0;
        Atotal(ubasis_function_number+1,:)=0;
        Atotal(ubasis_function_number+1,ubasis_function_number+1)=1;
        right(ubasis_function_number+1,1)=node_coordinate(1,6);
 

%%%%%%%%%%%%%%%%%% solution %%%%%%%%%%%%%%
%basisfunctionweight=Atotal\right;
pre=3;post=3;
cycle=1;
smooth=1;
grids=7;
maxit=1e6;
tol=1e-7;
%BT=B.';
IA=speye(size(A));
IB=speye(pbasis_function_number,pbasis_function_number);
D=diag(A);
Dainv=sparse(diag(1./D));
alpha=80;
 matprecond=[IA,-alpha*Dainv*B;sparse(pbasis_function_number,ubasis_function_number),IB];
Ahat=Atotal*matprecond;

[x,res] = multigridV2(ubasis_function_number,pbasis_function_number,Ahat,right,pre,post,cycle,smooth,grids,maxit,tol);
x=matprecond*x;
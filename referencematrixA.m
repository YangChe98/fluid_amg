function [A11reference,A12reference,A22reference]=referencematrixA()

syms xi eta;

A11reference=zeros(9,9);
A12reference=zeros(9,9);
A22reference=zeros(9,9);


% reference matrix A 
for i=1:9
       [phii,phiix,phiiy]=Q2basisfunction(i,xi,eta);
    for j=1:9
     
        [phij,phijx,phijy]=Q2basisfunction(j,xi,eta);
        A11=int(matlabFunction(phiix.*phijx,'vars',[xi eta]),xi,-1,1);
        A1=int(matlabFunction(A11,'var',[eta]),eta,-1,1);
        A11reference(i,j)=double(A1);
%         
      

        A22=int(matlabFunction(phiiy.*phijy,'vars',[xi eta]),xi,-1,1);
        A2=int(matlabFunction(A22,'var',[eta]),eta,-1,1);
        A22reference(i,j)=double(A2);
  
          A12=int(matlabFunction(phiix.*phijy,'vars',[xi eta]),xi,-1,1);
        A3=int(matlabFunction(A12,'var',[eta]),eta,-1,1);
        A12reference(i,j)=double(A3);


%         C11=int(matlabFunction(phix(2,i).*phix(2,j),'vars',[x y]),x,-1,1);
%        C1=int(matlabFunction(C11,'var',[y]),y,-1,1);
%          C11reference(i,j)=double(C1);
%          
%         C22=int(matlabFunction(phiy(2,i).*phiy(2,j),'vars',[x y]),x,-1,1);
%         C2=int(matlabFunction(C22,'var',[y]),y,-1,1);
%         C22reference(i,j)=double(C2);

% 
%         C12=int(matlabFunction(phix(2,i).*phiy(2,j),'vars',[x y]),x,-1,1);
%         C3=int(matlabFunction(C12,'var',[y]),y,-1,1);
%         C12reference(i,j)=double(C3);

    end
end

save A11reference.mat A11reference;
save A12reference.mat A12reference;
save A22reference.mat A22reference;
% save C11reference.mat C11reference;
% save C12reference.mat C12reference;
% save C22reference.mat C22reference;
end
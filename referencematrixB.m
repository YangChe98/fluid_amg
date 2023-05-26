function [B11reference,B12reference]=referencematrixB()

syms xi eta;

B11reference=zeros(9,4);
B12reference=zeros(9,4);
B21reference=zeros(9,4);
B22reference=zeros(9,4);
for i=1:9
    [phi,phix,phiy]=Q2basisfunction(i,xi,eta);
   
    
    for j=1:4
    
    
    [p,px,py]=Q1basisfunction(j,xi,eta);
    
        B11=int(matlabFunction(phix.*p,'vars',[xi eta]),xi,-1,1);
        B1=int(matlabFunction(B11,'var',[eta]),eta,-1,1);
        B11reference(i,j)=-double(B1);

         B12=int(matlabFunction(phiy.*p,'vars',[xi eta]),xi,-1,1);
        B2=int(matlabFunction(B12,'var',[eta]),eta,-1,1);
        B12reference(i,j)=-double(B2);


end
end

save B11reference.mat B11reference;
save B12reference.mat B12reference;

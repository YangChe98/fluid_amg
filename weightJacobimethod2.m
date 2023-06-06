function [x,res]=weightJacobimethod2(A,b,x0,weight,iter)

x=x0;
D=diag(A);
Dinv=diag(1./D);
Dinv=sparse(Dinv);
x_ana=A\b;
tol=1e-7;
norm_x_ana=norm(x_ana);
x_old=x;

for i=1:iter
    x=x+weight*Dinv*(b-A*x);
     abserror=norm(x-x_ana);
     res(i)=norm(x-x_ana)/norm_x_ana;
    if abserror<tol
        break;
    end
    x_old=x;
end
end


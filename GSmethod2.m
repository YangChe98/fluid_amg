function [x,res]=GSmethod2(A,b,x0,iter)

x=x0;
D=spdiags(diag(A));
U=triu(A,1);
L=tril(A,-1);
DLf=(D+L)\b;
DLU=(D+L)\U;
x_ana=A\b;
tol=1e-7;
norm_x_ana=norm(x_ana);
x_old=x;



for i=1:iter
    x=DLf-DLU*x;
     abserror=norm(x-x_ana)
     res(i)=norm(x-x_ana)/norm_x_ana;
    if abserror<tol
        break;
    end
end


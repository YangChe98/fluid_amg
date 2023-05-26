function rh=restrict(A,b,x,Amatrix,n)
r=b-A*x;
rh=Amatrix.R{n}*r;
end
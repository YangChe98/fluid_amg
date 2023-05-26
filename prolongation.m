function x = prolongation(x,Amatrix,n,eh)
x = x+Amatrix.P{n}*eh;
end

function x=smoothing(A,b,x,pre,smooth);
weight=0.1;
 if smooth==1
             x=weightJacobimethod(A,b,x,weight,pre);
         elseif smooth == 2 
               x=GSmethod(A,b,x,pre);
 end
end
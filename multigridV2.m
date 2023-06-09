function [x,res] = multigridV2(ubasis_function_number,pbasis_function_number,A,b,pre,post,cycle,smooth,grids,maxit,tol,x0)

%    A = A matrix (n x n)
%    b = Right hand side vector (n x 1)
%    pre = Number of presmoothing iterations
%    post = Number of postsmoothing iterations
%    cycle = Type of multigrid cycle (1=V-cycle, 2=W-cycle, 3=F-cycle)
%    smooth = Smoother type (1=Jacobi, 2=Gauss-Seidel)
%    grids = Max grids in cycle, used for grids level during recursion
%    maxit = Max iterations of solver
%    tol = Tolerance of solver
% OPTIONAL:
%    x0 = Initial Guess
% OUTPUTS:
%    x = Solution vector
%    res = Residual vector
    if  nargin<12
        x = b*0;
    else
        x = x0;
    end
    x_ana=A\b;
    norm_x_ana=norm(x_ana);
    Amatrix.Ac{1}=A;
Amatrix.A1{1}=A(1:ubasis_function_number,1:ubasis_function_number);
Amatrix.A2{1}=A(ubasis_function_number+1:ubasis_function_number+pbasis_function_number,ubasis_function_number+1:ubasis_function_number+pbasis_function_number);
     res = zeros(maxit,1);
     x_old=x;
    % weight=1; % weighted Jacobi method
     for i=1:grids
        [Amatrix.P1{i}, ~]=Beckcoarsen(Amatrix.A1{i});
         [Amatrix.P2{i}, ~]=Beckcoarsen(Amatrix.A2{i});
         [n1,m1]=size(Amatrix.P1{i});
          [n2,m2]=size(Amatrix.P2{i});
         Amatrix.P{i}=[Amatrix.P1{i},sparse(n1,m2);sparse(n2,m1),Amatrix.P2{i}];
         Amatrix.R{i}=Amatrix.P{i}.';
         Amatrix.A1{i+1}=Amatrix.P1{i}.'*Amatrix.A1{i}*Amatrix.P1{i};
          Amatrix.A2{i+1}=Amatrix.P2{i}.'*Amatrix.A2{i}*Amatrix.P2{i};
        Amatrix.Ac{i+1}=Amatrix.P{i}.'*Amatrix.Ac{i}*Amatrix.P{i};
    end
     res = zeros(maxit,1);
    % weight=1; % weighted Jacobi method
     for ii=1:maxit
         % presmoothing
%          if smooth==1
%              x=weightJacobimethod(A,b,x,weight,pre);
%          elseif smooth == 2 
%                x=GSmethod(A,b,x0,pre);
%          end
        Amatrix.f{1}=b;
        for k=1:grids
            if k==1
                Amatrix.x_tmp{k}=x;
            else
                Amatrix.x_tmp{k}=0* Amatrix.f{k};
            end
            if k==grids
%                 for j=1:pre^2
%                     Amatrix.x_tmp{k}=smoothing(Amatrix.Ac{k},Amatrix.f{k},Amatrix.x_tmp{k},pre,smooth);
%                 end
%                 Amatrix.rh{k}=Amatrix.x_tmp{k};
              Amatrix.rh{k}=Amatrix.Ac{k}\Amatrix.f{k};
            else

                for j=1:pre
                    Amatrix.x_tmp{k}=smoothing(Amatrix.Ac{k},Amatrix.f{k},Amatrix.x_tmp{k},pre,smooth);
                end
                Amatrix.rh{k}=restrict(Amatrix.Ac{k},Amatrix.f{k},Amatrix.x_tmp{k},Amatrix,k);
            end
            Amatrix.f{k+1}=Amatrix.rh{k};

        end


         for k=grids-1:-1:1
         
                Amatrix.x_tmp{k}=prolongation(Amatrix.x_tmp{k},Amatrix,k,Amatrix.x_tmp{k+1});
            
            for j=1:post
                    Amatrix.x_tmp{k}=smoothing(Amatrix.Ac{k},Amatrix.f{k},Amatrix.x_tmp{k},pre,smooth);
             end
        end
    x=Amatrix.x_tmp{1};
    abserror=norm(x-x_ana);
     res(ii)=norm(x-x_ana)/norm_x_ana;
    if abserror<tol
        break;
    end
    x_old=x;

     end
   end

         
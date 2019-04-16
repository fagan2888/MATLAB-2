%Joshua Enxing
%Tufts University
%MA150

%Conjugate Gradient Method
%Solves the equation Ax = b
%Inputs A,b,initial guess x_0, max number of iters, error threshold
%Returns approximate solution vector x, number of iterations to reach
%threshold
function [x,iters] = conjugate_gradient(A,b,x_0,max_iters,tol)

x = x_0;
r = b-A*x;
p = r;
iters = 0;

%stop when max number of iterations have been reached
for i=1:max_iters
    
    %stop when below tolerance
    while norm(r) >= tol
        
        %update alpha
        alpha = (r'*r)/(p'*A*p);
        
        %update x
        x = x + alpha*p;
        temp = r;
        
        %update r
        r = r - alpha*A*p;
        
        %update beta
        beta = (r'*r)/(temp'*temp);
        
        %update p
        p = r + beta*p;
        
        %add 1 to number of iterations
        iters = iters + 1;
    end
end

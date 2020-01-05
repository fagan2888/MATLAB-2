%Joshua Enxing
%Tufts University
%MA150

%Regular SGD
%Inputs grad_f vector of gradient functions, m number of terms in the sum
%for f, x0 initial guess, constant s for stepsize diminishing, max number of (outer) iterations
%max_iters, real solution sol, and tolerance tol
function [x,err] = sgd(grad_f,m,x0,s,max_iters,sol,tol,seed,R)

err = zeros(max_iters,1);
x = x0;
disp("Iteration || Norm of x_k - x^*");

%fix seed
rng(seed);

%For max_iters outer iterations, we do m*max_iters total iterations
for i=1:max_iters
    
    %Find absolute error and store in vector
    v = norm(x - sol);
    err(i) = v;
    disp(i + " || " + v);
    
    %if error is above tol, continue, else break
    if v > tol
        
        %inner loop
        for j=1:m
            grad_j = grad_f{randi(m)};
            
            %use alpha = 1/((i)^s)
            x = x - (R/((i)^s))*grad_j(x);
        end
    else
        err = err(1:i);
        break;
    end
end
end
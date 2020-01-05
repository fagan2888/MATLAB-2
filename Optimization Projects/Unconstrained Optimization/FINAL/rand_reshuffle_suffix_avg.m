%Joshua Enxing
%Tufts University
%MA150

%SGD with random reshuffling and q-suffix averaging
%Inputs grad_f vector of gradient functions, m number of terms in the sum
%for f, x0 initial guess, alpha stepsize, max number of (outer) iterations
%max_iters, real solution sol, and tolerance tol
function [x,err,iters] = rand_reshuffle_suffix_avg(grad_f,m,x0,s,q,max_iters,sol,tol,seed,R)

err = zeros(max_iters*q,1);
iters = zeros(length(err),1);
x = x0;
disp("Iteration || Norm of x_k - x^*");

%fix seed
rng(seed);

l = 1;
x_temp = x;
v = norm(x_temp - sol);
x_avg = zeros(length(x),max_iters);
%For max_iters outer iterations, we do m*max_iters total iterations
for i=1:max_iters
    
    %permute the elements 1,...,m
    p = randperm(m);
    
    %when qk is an integer
    if mod(q*i,1) == 0
        
        %Find absolute error and store in vector
        v = norm(x_temp - sol);
        err(l) = v;
        iters(l) = i;
        disp(i + " || " + v);
        l = l+1;
    end
    
    if v > tol
        %inner loop
        alpha = (R/((i+1)^s));
        for j=1:m
            grad_j = grad_f{p(j)};
            
            %use alpha_i = 1/((i+1)^s)
            x = x - alpha*grad_j(x);
        end
        
        %filling vector of xbar_{1,i}
        if i == 1
            x_avg(:,i) = x;
        else
            x_avg(:,i) = ((i-1)*x_avg(:,i-1) + x)/i;
        end
        
        %when qk is an integer
        if mod(q*i,1) == 0
            x_temp = (x_avg(:,i) - q*x_avg(:,round((1-q)*i)))/(1-q);
        end
    else
        err = err(1:l-1); 
        iters = iters(1:l-1);
        break;
    end
end
err = err(1:l-1);
iters = iters(1:l-1);
x = x_temp;
end
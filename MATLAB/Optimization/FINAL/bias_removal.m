%Joshua Enxing
%Tufts University
%MA150

%SGD with random reshuffling and q-suffix averaging with bias removal
%Inputs grad_f vector of gradient functions, m number of terms in the sum
%for f, x0 initial guess, constant s for stepsize diminishing, value q, max number of (outer) iterations
%max_iters, real solution sol, and tolerance tol, initial paramete mu_0,
%initial parameter H_0
function [x,err,iters] = bias_removal(grad_f,hess_f,m,x0,s,q,max_iters,sol,tol,mu_0,H_0,seed,R)

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
alpha_avg = zeros(max_iters,1);

mu = mu_0;
H = H_0;
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
            hess_j = hess_f{p(j)};
            
            if i == max_iters-1
                mu = mu + hess_j(x)*grad_j(x);
                H = H + hess_j(x);
            end
            
            %use alpha_i = 1/((i+1)^s)
            x = x - alpha*grad_j(x);     
        end
        
        %filling vector of xbar_{1,i}
        if i == 1
            x_avg(:,i) = x;
            alpha_avg(i) = R/(2^s);
        else
            x_avg(:,i) = ((i-1)*x_avg(:,i-1) + x)/i;
            alpha_avg(i) = ((i-1)*alpha_avg(i-1) + alpha)/i;
        end
        
        %when qk is an integer
        if mod(q*i,1) == 0
            x_temp = (x_avg(:,i) - q*x_avg(:,round((1-q)*i)))/(1-q);
        end
        
        if i == max_iters
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

alpha_temp = (alpha_avg(max_iters)-q*alpha_avg(round((1-q)*max_iters)))/(1-q);
%bqk = -(alpha_temp*mu)/H;
bqk = -H\(alpha_temp*mu);
x = x_temp - bqk;
end
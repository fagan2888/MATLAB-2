%Projected newton method

%Inputs: function to minimize f, gradient of f grad_f, hessian of f hess_f,
%l vector of lower bound constraints, u vector of upper bound constraints, 
%max iterations max_iters, tolerance tol, alpha value alpha, beta value beta
function [y,x] = proj_newt(x_0,f,grad_f,hess_f,u,l,max_iters,tol,alpha,beta)

x = x_0;

iters = 0;
x_last = x - ones(size(x));
err = norm(x-x_last);

y = zeros(max_iters,2);
y(1,:) = x';

df = grad_f(x(1),x(2)); 
eps = min((u-l)/2);
while iters < max_iters && err > tol
    %Find reduced hessian
    df = grad_f(x(1),x(2)); 
    R = reduced_hess(x,eps,u,l,hess_f); 
    
    %Check if reduced hessian is SPD
    [~,p] = chol(R);
    if p~=0 || norm(R'-R) > 0
        disp("Error: reduced Hessian is not SPD");
        break;
    end
    
    %search direction
    d = -R\df; 
    
    %find stepsize
    lambda = find_lambda(alpha,df,beta,f,x,d,u,l);
    
    %calculate eps in projected newton
    eps = min(norm(x - project(x+d,u,l)),min(u-l)/2);
    err = eps;
    
    %update iterate
    x = project(x + lambda*d,u,l);
    
    %keep track of iterates
    y(iters+2,:) = x';
    iters = iters+1;
end

y = y(1:iters,:);
    

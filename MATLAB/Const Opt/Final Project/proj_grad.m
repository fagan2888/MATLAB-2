%Projected gradient method with armijo rule 

%Inputs: initial guess z_0, function f, gradient df, upper and lower bound
%constraints u and l, alpha_0 c and p for armijo, max number iterations
%max_iters, and tolerance tol
function [y,z] = proj_grad(z_0, f, df, u, l, alpha_0, p, c, max_iters, tol)

z = z_0;
alpha = alpha_0;
y = zeros(1,2);
y(1,:) = z';
err = 1;
%break condition
iter = 1;
while iter < max_iters && err > tol
    %descent direction
    p_k = -1*df(z(1),z(2));
   
    d = -p_k;
    
    %armijo rule
    fxk = f(z(1),z(2));
    while f(z(1)+alpha*p_k(1),z(2)+alpha*p_k(2)) > fxk + c*alpha*p_k'*d
        alpha = p*alpha;
    end
    
    err = norm(z - project(z+p_k,u,l));
    %update guess
    z = z + alpha*p_k;
    
    %project if necessary
    z = project(z,u,l);
    y(iter+1,:) = z';
    iter = iter+1;
end

y = y(1:iter,:);


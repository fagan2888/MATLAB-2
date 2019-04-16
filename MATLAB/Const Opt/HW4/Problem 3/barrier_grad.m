%Barrier method for gradient descent, alpha determined by armijo rule with
%backtraking

%Inputs: x initial guess, parameter mu, tolerance tol, max number of
%iterations max_iters, parameter alpha, parameter rho, parameter c

%Suggested inputs: x [1 1]', mu value 1e-m, where m is a
%nonnegative integer, max_iters 100, alpha = rho = c = 0.5

%Outputs: convergence history iters and final answer x
function [iters,x] = barrier_grad(x,mu,tol,max_iters,alpha_0,rho,c)

%Objective function
f = @(x,y) 0.5*(x^(2)+y^(2));

%Gradient of f
gradf = @(x,y) [x;y];

%Penalty function
phi = @(x,y) -mu*log(x^(2)+y^(2) - 1);

%Gradient of penalty function
grad_phi = @(x,y) -2*mu*[x/(x^(2)+y^(2)-1); y/(x^(2)+y^(2)-1)];

%Barrier function
B = @(x,y) f(x,y) + phi(x,y); 

%Gradient of barrier function
gradB = @(x,y) gradf(x,y) + grad_phi(x,y);

% Keep track of iterate history:
iters = zeros(max_iters,2);
err = zeros(max_iters+1,1);
lim_err = zeros(max_iters,1);
err(1) = norm(x);

x_last = [1000 1000]';
iter = 1;
% Apply unconstrained Gradient Descent to penalty function B(x,y):
while iter < max_iters && norm(x_last - x) > tol
    %descent direction:
    pk = -1*gradB(x(1),x(2));
    
    alpha = alpha_0;
    %armijo rule with backtracking
    while(B(x(1)+alpha*pk(1),x(2)+alpha*pk(2))>B(x(1),x(2))+c*alpha*pk'*gradB(x(1),x(2)))|| (x(1)+alpha*pk(1) + x(2)+alpha*pk(2) < 1)
        alpha = rho*alpha;
    end
    
    %update iterate list
    iters(iter,:) = x;
    
    %update x
    x_last = x;
    x = x + alpha*pk;
    
    iter = iter+1;
end
iters = iters(1:iter-1,:);


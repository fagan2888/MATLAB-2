%Newton's method with merit function

%Inputs: x, lambda
function [x,hist] = newton_merit(x, lambda, alpha, rho, c, max_iters, tol)

%f,g,and M
f = @(x,y) exp(3*x + 4*y);
g = @(x,y) x^2 + y^2 - 1;
M = @(x,y) f(x,y) + 100*norm(g(x,y))^2;

%grad f, grad g, grad M and hess f, hess g
f_x = @(x,y) exp(3*x + 4*y)*[3;4];
g_x = @(x,y) 2*[x;y];
M_x  = @(x,y)  f_x(x,y) + 200* g(x,y)* g_x(x,y);

f_xx = @(x,y) exp(3*x + 4*y)* [9, 12; 12, 16];
g_xx =  @(x,y)[2, 0; 0, 2];

%grad lagrangian and hess lagrangian
L_x = @(x,y, lambda) f_x(x,y) - lambda*g_x(x,y);
L_xx = @(x,y, lambda) f_xx(x,y) - lambda* g_xx(x,y);

%newton matrices
H = @(x,y, lambda) [L_xx(x,y, lambda), -g_x(x,y);  -(g_x(x,y))', 0];
b = @(x,y, lambda) [-L_x(x,y, lambda); g(x,y)];

x1 = x(1);
x2 = x(2);
Lambda = lambda;
Merit = M(x(1),x(2));

%newton's method
for k = 1: max_iters
    x_prev = x;
    lambda_prev = lambda;
    
    %descent direction
    p = H(x(1), x(2), lambda) \ b(x(1), x(2), lambda);
    
    %line search
    while M(x(1) + alpha*p(1), x(2) + alpha*p(2)) > M(x(1), x(2)) + c*alpha*(M_x(x(1), x(2)))'*p(1:2)
        alpha = rho*alpha;
    end
    
    %update
    x = x + alpha* p(1:2);
    lambda = lambda + alpha* p(3);
    merit = M(x(1),x(2));
    
    %convergence history
    x1 = [x1; x(1)];
    x2 = [x2; x(2)];
    Lambda = [Lambda; lambda];
    Merit = [Merit; merit];
    
    %stopping criterion
    if norm (x - x_prev) < tol && norm (lambda - lambda_prev) < tol
        break
    end
    
end

hist = [x1 x2 Lambda Merit];
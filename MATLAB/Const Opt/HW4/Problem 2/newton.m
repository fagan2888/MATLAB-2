%Newton's method without merit function

%Inputs: x, lambda, max_iters, tolerance

%Outputs: solution x, 
function  [x,hist] = newton(x,lambda,max_iters,tol)

%f and g
f = @(x,y) exp(3*x + 4*y);
g = @(x,y) x^2 + y^2 - 1;

%grad f and grad g
f_x = @(x,y) exp(3*x + 4*y)*[3;4];
g_x = @(x,y) 2*[x;y];

%hessian f and hessian g 
f_xx = @(x,y) exp(3*x + 4*y)* [9, 12; 12, 16];
g_xx =  @(x,y)[2, 0; 0, 2];

%lagrangian gradient and hessian
L_x = @(x,y, lambda) f_x(x,y) - lambda*g_x(x,y);
L_xx = @(x,y, lambda) f_xx(x,y) - lambda* g_xx(x,y);

%matrices for newtons method
H = @(x,y, lambda) [L_xx(x,y, lambda), -g_x(x,y);  -(g_x(x,y))', 0];
b = @(x,y, lambda) [-L_x(x,y, lambda); g(x,y)];

x1 = x(1);
x2 = x(2);
Lambda = lambda;

%newton's method
for k = 1: max_iters
    x_prev = x;
    lambda_prev = lambda;
    
    p = H(x(1), x(2), lambda) \ b(x(1), x(2), lambda);
    x = x + p(1:2);
    lambda = lambda + p(3);
    
    %convergence history
    x1 = [x1; x(1)];
    x2 = [x2; x(2)];
    Lambda = [Lambda; lambda];
    
    %if close enough, stop
    if norm (x - x_prev) < tol && norm (lambda - lambda_prev) < tol
        break
    end
    
end

%convergence history
hist = [x1 x2 Lambda];

    
    
    
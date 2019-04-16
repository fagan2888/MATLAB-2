%function, gradient, and hessian
 f = @(x,y) (x^2-1)^2 + (x^2*y-x-1)^2;
 df = @(x,y)[(4*x*(x^2-1)+(4*x*y-2)*(x^2*y-x-1)) 2*x^2*(x^2*y-x-1)]';
 hess_f = @(x,y)[(12*x^2-2+12*x^2*y^2-12*x*y-4*y) (8*x^3*y-6*x^2-4*x);(8*x^3*y-6*x^2-4*x) 2*x^4];

%set values
tol = 1E-10;
alpha_newt = 1E-4;
beta = 0.1;
max_iters = 100;

alpha_grad = 0.5;
p = 0.5;
c = 0.5;

%set initial guess and bound constraints
z_0 = [0 0]';
l = [-1 1]'; 
u = [0 2]';

[y1,z1] = proj_newt(z_0,f,df,hess_f,u,l,max_iters,tol,alpha_newt,beta);

[y2,z2] = proj_grad(z_0, f, df, u, l, alpha_grad, p, c, max_iters, tol);
%function, gradient, and hessian
% f = @(x,y) x^2 - y^2;
% df = @(x,y)[2*x -2*y]';
% hess_f = @(x,y)[2 0;0 -2];


%set values
tol = 1E-6;
alpha_newt = 1E-4;
beta = 0.9;
max_iters = 100;

alpha_grad = 0.5;
p = 0.5;
c = 0.5;

%set initial guess and bound constraints
% z_0 = [0.2 0.2]';
% l = [-1 -1]'; 
% u = [1 1]';
% 
% [y1,z1] = proj_newt(z_0,f,df,hess_f,u,l,max_iters,tol,alpha_newt,beta);
% 
% [y2,z2] = proj_grad(z_0, f, df, u, l, alpha_grad, p, c, max_iters, tol);

%set initial guess and bound constraints
z_0 = [0.2 0.2]';
l = [0 0]'; 
u = [2 2]';

[y1,z1] = proj_newt(z_0,f,df,hess_f,u,l,max_iters,tol,alpha_newt,beta);

[y2,z2] = proj_grad(z_0, f, df, u, l, alpha_grad, p, c, max_iters, tol);
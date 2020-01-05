%Joshua Enxing
%Tufts University
%MA150

%n-dimensional Newton's Method with armijo rule. 

%Takes as inputs:

%symbolic function f (must first type syms x1 x2 ... xn into terminal),
%vector of variables vars (e.g. [x y]),
%initial guess vector x_0, 
%intial alpha alpha_0
%rho for the armijo rule,
%c for the armijo rule
%max number of iterations max_iters,
%tolerance tol of the norm of the gradient

function [z,iters] = newtons_method_ndim_armijo(f,vars,x_0,alpha_0,rho,c,max_iters,tol)

syms x y;
z = x_0;
iters = 0;
grad = gradient(f);
hess = hessian(f);

g = zeros(length(z),1);
h_mat = zeros(length(g));

alpha = alpha_0;
for i=1:max_iters
    
    for j=1:length(g)
        g(j) = subs(grad(j),vars,z');
    end
    
    disp(iters + " " + norm(g));
    %stop if norm of gradient is smaller than tol
    if norm(g) < tol
        break;
    end
    
    for k=1:length(g)
        for l=1:length(g)
            h_mat(k,l) = subs(hess(k,l),vars,z');
        end
    end
    
    p = -h_mat\g;
    
    fxk = subs(f,vars,z');
    while subs(f,vars,(z+alpha*p)') > fxk + c*alpha*p*g'
        alpha = rho*alpha;
    end
    
    z = z + alpha*p;
    iters = iters + 1;
end



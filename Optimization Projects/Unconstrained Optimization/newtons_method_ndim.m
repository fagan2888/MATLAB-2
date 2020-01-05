%Joshua Enxing
%Tufts University
%MA150

%n-dimensional Newton's Method. Takes symbolic function f,
%must first type syms x1 x2 ... xn into terminal,
%vector of variables vars (e.g. [x y])
%initial guess vector x_0, max number of iterations max_iters,
%tolerance tol of the norm of the gradient

function [z,iters] = newtons_method_ndim(f,vars,x_0,max_iters,tol)

syms x y;
z = x_0;
iters = 0;
grad = gradient(f);
hess = hessian(f);

g = zeros(length(z),1);
h_mat = zeros(length(g));
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
    z = z + p;
    iters = iters + 1;
end



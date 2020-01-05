%Joshua Enxing
%Tufts University
%MA150

%n-dimensional Trust-Region Method with Dogleg. Takes symbolic function f,
%must first type syms x1 x2 ... xn into terminal,
%vector of variables vars (e.g. [x y])
%initial guess vector x_0, 
%initial guess radius r_0,
%nu in dogleg method,
%max number of iterations max_iters,
%tolerance tol of the norm of the gradient

function [z,iters] = trust_region_dogleg(f,vars,x_0,r_0,nu,max_iters,tol)

disp("x_0: " + x_0(1) + ", y_0: " + x_0(2) + ", r_0: " + r_0 + ", eta: " + nu);  
disp("Iteration || Norm of gradient");
syms x y;
z = x_0;
r = r_0;
iters = 0;

grad = gradient(f);
hess = hessian(f);

fxk = subs(f,vars,z');
g = evaluate_gradient(grad,vars,z);
h_mat = evaluate_hessian(hess,vars,z);

pkqn = -h_mat\g;
for i=1:max_iters
    
    pk = determine_pk(g,h_mat,pkqn,r);
    fxk_plus_pk = subs(f,vars,(z+pk)');
    mkpk = fxk + (g')*pk + 0.5*(pk')*h_mat*pk;
    
    rho = (fxk - fxk_plus_pk)/(fxk - mkpk);
    
    if rho >= nu
        z = z + pk;
        fxk = subs(f,vars,z');
        g = evaluate_gradient(grad,vars,z);
        
        h_mat = evaluate_hessian(hess,vars,z);
        pkqn = -h_mat\g;
    end
    
    disp(iters + " || " + norm(g));
    %stop if norm of gradient is smaller than tol
    if norm(g) < tol
        break;
    end
    
    if rho > 0.75
        r = min(2*r,r_0);
    end
    if rho < 0.25
        r = 0.25*r;
    end
    
    iters = iters + 1;
end



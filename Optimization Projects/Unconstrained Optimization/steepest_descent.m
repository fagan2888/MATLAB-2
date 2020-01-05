%Joshua Enxing
%Tufts University
%MA150

%n-dimensional steepest descent with armijo rule. Takes symbolic function f,
%must first type syms x1 x2 ... xn into terminal,
%vector of variables vars (e.g. [x y])
%initial guess vector z_0, alpha_0 initial alpha for armijo rule,
%p for armijo rule, c for armijo rule, max_iters max # of iterations,
%and tolerance tol
function z = steepest_descent_armijo(f, vars, z_0, alpha_0, p, c, max_iters, tol)

iters = 0;
z = z_0;
alpha = alpha_0;
syms x y
df = gradient(f,[x y]);
for i=1:max_iters
    p_k = [double(-1*subs(df(1),vars,z')) double(-1*subs(df(2),vars,z'))];
    d = -p_k;
    
    if norm(p_k) < tol
        break;
    end
    
    while double(subs(f,vars,z'+alpha*p_k)) > double(subs(f,vars,z'+alpha*p_k)) + c*alpha*p_k*d'
        alpha = p*alpha;
    end
    
    z = z + alpha*p_k';
end


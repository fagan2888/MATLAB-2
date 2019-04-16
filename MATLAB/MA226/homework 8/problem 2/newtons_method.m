function p = newtons_method(p0,max_iters,tol)

p = p0;
f = f_eval(p(1),p(2),p(3),p(4));
disp("Iteration || Norm residual");

for i=1:max_iters
    disp(i + " || " + norm(f));
    if norm(f) > tol
        J = jacobian(p(1),p(2),p(3),p(4));
        delta = J\f;
        p = p - delta;
        f = f_eval(p(1),p(2),p(3),p(4));
    else 
        break;
    end
end
function [x,i] = steep_desc(A,b,c,x0,tol,max_iters)

x = x0;
%disp("Iteration || Norm of Gradient");
grad = grad_quad(A,b,x);
for i=1:max_iters
    %disp(i - 1 + " || " + norm(grad));
    if norm(grad) > tol
        s = -grad;
        alpha = compute_alpha(A,b,s,x);
        x = x + alpha*s;
        grad = grad_quad(A,b,x);
    else
        break;
    end
end
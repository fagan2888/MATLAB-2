function [x,err] = gradient_descent(A,b,m,x0,alpha,max_iters,sol,tol)

err = zeros(max_iters,1);
x = x0;
disp("Iteration || Norm of x_k - x^*");
[grad,grad_full] = grad_eval(A,b,x,m);
for i=1:max_iters
    v = norm(x-sol);
    err(i) = v;
    disp(i-1 + " || " + v);
    if v > tol
        x = x - alpha*grad_full;
        [grad,grad_full] = grad_eval(A,b,x,m);
    else
        err = err(1:i);
        break;
    end
end


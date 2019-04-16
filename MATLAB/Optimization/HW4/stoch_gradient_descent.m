function [x,err] = stoch_gradient_descent(A,b,m,x0,alpha,max_iters,sol,tol)

err = zeros(max_iters,1);
x = x0;
disp("Iteration || Norm of x_k - x^*");
for i=1:max_iters
    v = norm(x - sol);
    err(i) = v;
    disp(i-1 + " || " + v);
    if v > tol
        j = randi(m);
        grad = (1/10)*(A{j}*x -b{j});
        x = x - alpha*grad;
    else
        err = err(1:i);
        break;
    end
end
function [x,err] = reshuffling_gradient_descent(A,b,m,x0,alpha,max_iters,sol,tol)

err = zeros(max_iters,1);
x = x0;
disp("Iteration || Norm of x_k - x^*");
for i=1:max_iters
    if mod(i,m) == 1
        p = randperm(10);
    end
    v = norm(x - sol);
    err(i) = v;
    disp(i-1 + " || " + v);
    
    
    if mod(i,m) == 0
        j = p(10);
    else
        j = p(mod(i,m));
    end
    
    if v > tol
        grad = (1/10)*(A{j}*x -b{j});
        x = x - alpha*grad;
    else
        err = err(1:i);
        break;
    end
end
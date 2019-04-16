function [r,rel,e] = run_weighted_jacobi_tridiag(A,x_guess,x_tru,b_vect,w_weight,max_iters)

w = w_weight;
x = x_guess;
b = b_vect;
r = zeros(max_iters,1);
rel = zeros(max_iters-1,1);
e = zeros(max_iters,1);

for l=1:max_iters
    x = weighted_jacobi_tridiag(A,x,b,w);
    r(l) = norm(b - A*x);
    if l>1
        rel(l) = r(l)/r(l-1);
    end
    e(l) = norm(x - x_tru);
end

close all;

figure
plot(rel);
xlabel('Number of Iterations');
ylabel('Norm of Relative Residual Vector');
title('Weighted Jacobi Norm of Relative Residuals');

figure
plot(r);
xlabel('Number of Iterations');
ylabel('Norm of Residual Vector');
title('Weighted Jacobi Norm of Residuals');

figure
plot(e);
xlabel('Number of Iterations');
ylabel('Norm of Error Vector');
title('Weighted Jacobi Norm of Error');
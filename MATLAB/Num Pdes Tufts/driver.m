x_guess = rand(129,1);
b = zeros(129,1);
A = hw_5_make_matrix_A(129);

[x,error,res,rel_res] = run_vcycle_tridiag(x_guess,b,A,b,1,1,200);

figure
plot(rel_res);
xlabel('Number of Iterations');
ylabel('Norm of Relative Residual Vector');
title('V-cycle Norm of Relative Residuals, k = 200, n = 129');

figure
plot(res);
xlabel('Number of Iterations');
ylabel('Norm of Residual Vector');
title('V-cycle Norm of Residuals, k = 200, n = 129');

figure
plot(error);
xlabel('Number of Iterations');
ylabel('Norm of Error Vector');
title('V-cycle Norm of Error, k = 200, n = 129');
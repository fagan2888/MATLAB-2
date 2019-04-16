%One iteration of two-grid scheme for tridiagonal matrix A
function x = two_grid_tridiag(A, x_guess, b_vec, nu_1)

n = length(x_guess);
x = x_guess;

%Makes the interpolation matrix P
P = make_matrix_P(n);

%Relax nu_1 times with gauss-seidel from hw5
for j=1:nu_1
    x = gauss_seidel_tridiag(A,x,b_vec);
end

%Restrict residual
b = P'*(b_vec-A*x);

%Correct
e = (P'*A*P)\b;
x = x + P*e;



    
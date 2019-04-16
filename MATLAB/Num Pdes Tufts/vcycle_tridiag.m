%Does one iteration of vcycle
function x = vcycle_tridiag(x_h,b_h,A_h,nu_1,nu_2)

x = x_h;
b = b_h;
A = A_h;
n = length(x);

%Relax nu_1 times with gauss-seidel from hw5
for i=1:nu_1
    x = gauss_seidel_tridiag(A,x,b);
end

%make interpolation matrix P
P = make_matrix_P(n);

%Restrict A
A_2h = (P')*A*P;

%Restrict residual
res = (P')*(b - A*x);

%If on next-to-coarsest grid, restrict and direct-solve, else
%recursively call vcycle solving Ae = r on the coarser grid
if n == 5
    e = A_2h\res;
else 
    temp = zeros(length(res),1);
    e = vcycle_tridiag(temp,res,A_2h,nu_1,nu_2);
end

%Correct
x = x + P*e;

%Relax nu_2 times with gauss-seidel from hw5
for i=1:nu_2
    x = gauss_seidel_tridiag(A,x,b);
end

    
    
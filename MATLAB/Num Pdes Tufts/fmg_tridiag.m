function x = fmg_tridiag(b_h,A_h,nu_1,nu_2)
disp('fmg');

A = A_h;
b = b_h;
n = length(b);
temp = zeros(n,1);
P = make_matrix_P(n);

if length(b) == 3
    x = A\b;
else
    A_2h = (P')*A*P;
    b = P'*b;
    x = fmg_tridiag(b,A_2h,nu_1,nu_2);
end

n = length(x);
P = make_matrix_P(2*(n-1)+1);
x = P*x;
b = P*b;
A = P*A*(P');

x = vcycle_tridiag(x,b,A,nu_1,nu_2);
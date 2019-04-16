%Makes the finite-element matrix for a partition of [0,1] of 
%num_unknowns + 2 points using piecewise linears for -u'' = f with 
%Dirichlet BCs
function A = hw_5_make_matrix_A(num_unknowns)

n = num_unknowns;
h = 1/(n+1);

%Make matrix A as described in lecture
diag = zeros(n,1);
offdiag = zeros(n-1,1);
for i=1:n-1
    diag(i) = 2/h;
    offdiag(i) = (-1)*(1/h);
end
diag(n) = 2/h;

A = gallery('tridiag',offdiag,diag,offdiag);

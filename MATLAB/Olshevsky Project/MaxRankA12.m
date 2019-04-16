%Find maximum rank of the upper-right corner of each symmetric 
%partition of the matrix into four submatrices.

function m = MaxRankA12(C)

[m,n] = size(C);
m = 0;

%Find rank of each of the n-1 submatrices of interest
for i=1:n-1
    
    %For each i, A is matrix of the first i rows and last n-i columns
    disp(C(1:i,i+1:n));
    
    %Matlab's rank algorithm 
    s = svd(C(1:i,i+1:n));
    disp(s);
    tol = max(size(C(1:i,i+1:n)))*eps(max(s));
    r = sum(s > tol);
    
    disp(r);
    
    %m is the maximum of the rank of A and the ranks of all other 
    %previous submatrices considered
    m = max(m,r);
end
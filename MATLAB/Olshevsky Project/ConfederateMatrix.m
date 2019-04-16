%Given a moment matrix M, we find matrix of orthogonal polynomials U
%corresponding to M. We then find the confederate matrix of U and 
%call it C.

function [C,U] = ConfederateMatrix(M)

[m,n] = size(M);
U = zeros(n);

%Find columns of U, the orthogonal polynomials related to the
%moment matrix M
for k=1:n
    v = zeros(k,1);
    v(k) = 1;
    
    MK = M(1:k,1:k);
    u = MK\v;
    U(1:k,k) = u;
end

C = zeros(n-1);

%Find (n-1)x(n-1) confederate matrix of the polynomials associated 
%with U
for i=1:n-1
    if i ~= n-1
        %N is (i+1)x(i+1) leading submatrix of U
        N = U(1:i+1,1:i+1);
        v = zeros(i+1,1);
        
        % v is x*u_(k-1)
        v(2:i+1) = U(1:i,i);
        x = N\v;
        C(1:i+1,i) = x;
    else
        %when i = n-1, N is ixi leading submatrix of U
        N = U(1:i+1,1:i+1);
        v = zeros(i+1,1);
        
        %Now use only first i entries of x, as opposed to first
        %i+1 entries as before
        v(2:i+1) = U(1:i,i);
        x = N\v;
        C(1:i,i) = x(1:i); 
    end
end


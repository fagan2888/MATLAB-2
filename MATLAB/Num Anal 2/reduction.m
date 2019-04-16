%Joshua Enxing
%University of Connecticut
%MA5511

%Reduces a real symmetric matrix to tridiagonal symmetric form
function [T,F] = reduction(S)

[m,n] = size(S);

%Create array to hold diagonal elements of reduced matrix
diag = zeros(n,1);

%Create array to hold subdiagonal elements of reduced matrix, also
%superdiagonal elements since symmetric
subdiag = zeros(n-1,1);

F = eye(n);

%Each iteration zeros out entries below the subdiagonal, and hence
%above the superdiagonal, since we are making similarity transformations
%that preserve the symmetry
for j=1:n-2
    diag(j) = S(j,j);
    s = norm(S(j+1:n,j));
    a = S(j+1,j);
    absa = abs(a);
    
    %Make sure we don't divide by something near 0
    if (s < 10^(-8))
        b = 0;
    else
        b = 1/(s*(s+absa));
    end
    
    %Make sure we don't divide by something near 0
    if (absa < 10^(-8))
        k = 0;
    else
        k = -s*(a/absa);
    end
    
    subdiag(j) = k;
    
    %Vector u
    u = S(j+1:n,j);
    u(1) = u(1)-k;
    
    %Matrix L
    L = S(j+1:n,j+1:n);
    
    %Vector p
    p = b*L*u;
    
    %Make q
    q = p -(b/2)*(p'*u)*u;
    
    S(j,j) = b;
    S(j+1,j) = a-k;
    
    %Fill in bottom right submatrix, what we want for next iteration
    S(j+1:n,j+1:n) = L - q*u'- u*q';
    S(j+1:n,j) = u;
    
    %Form the householder matrix
    T = eye(n);
    T(j+1:n,j+1:n) = T(j+1:n,j+1:n) - b*(u*u');
    
    %Compose all of the householder matrices so we can
    %use to find eigenvectors later
    F = F*T;
    
end
diag(n-1) = S(n-1,n-1);
diag(n) = S(n,n);
subdiag(n-1) = S(n,n-1);

T = eye(n);

%Fill in desired tridiagonal matrix
for i=1:n-1
    T(i,i) = diag(i);
    T(i+1,i) = subdiag(i);
    T(i,i+1) = subdiag(i);
end

T(n,n) = diag(n);
    
%Joshua Enxing
%University of Connecticut
%MA5511

%Does one iterate of the QR algorithm for a symmetric tridiagonal
%matrix A
function [Q,A] = qrIterateTridiag(A)

[n,m] = size(A);
H = eye(n);
    %The jth step forms a givens matrix to annihilate the entry in the
    %jth column and (j+1)st row
    for j=1:n-1
        %Find values for givens matrix
        if (((A(j,j))^2 + (A(j+1,j))^2) ~= 0)
            c = A(j,j)/sqrt(A(j,j)^2 + A(j+1,j)^2);
            s = A(j+1,j)/sqrt(A(j,j)^2 + A(j+1,j)^2);
        else
            c = 1;
            s = 0;
        end
        
        %Make givens matrix
        G = eye(2);
        G(1,1) = c;
        G(1,2) = s;
        G(2,1) = -s;
        G(2,2) = c;
        
        %Left multiply A by G
        A(j:j+1,1:n) = G*A(j:j+1,1:n);
        
        G(1,2) = -s;
        G(2,1) = s;
        
        %Right multiply A by G'
        A(1:n,j:j+1) = A(1:n,j:j+1)*G;
        H(1:n,j:j+1) = H(1:n,j:j+1)*G;
    end
% A is the result of the iteration, and Q is the corresponding Q, 
% so we can back-transform to get eigenvectors
Q = H;
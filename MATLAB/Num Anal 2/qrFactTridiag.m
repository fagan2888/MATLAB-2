function [Q,A] = qrIterateTridiag(A)

[n,m] = size(A);
H = eye(n);
    for j=1:n-1
        c = abs(A(j,j))/sqrt((A(j,j))^2 + (A(j+1,j))^2);
        s = abs(A(j+1,j))/sqrt((A(j,j))^2 + (A(j+1,j))^2);
        
        if (A(j,j)*A(j+1,j) ~= 0)
            k = (A(j,j)*A(j+1,j))/abs(A(j,j)*A(j+1,j));
        else
            k = 0;
        end
        
        G = eye(n);
        G(j,j) = c;
        G(j,j+1) = k*s;
        G(j+1,j) = -k*s;
        G(j+1,j+1) = c;
        
        %left mult
        A = G*A;
        
        G(j,j+1) = -k*s;
        G(j+1,j) = k*s;
        %right mult
        A = A*G;
        H = H*G;
    end
Q = H;
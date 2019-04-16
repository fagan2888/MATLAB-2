%Joshua Enxing
%University of Connecticut
%MA5511

%QR algorithm  on a real, symmetric, nonsingular matrix S
%with Rayleigh quotient shift, shifting by the 
%bottom right entry of the current iterate
function symmetricQRshift(S)

%Reduce S to symmetric tridiagonal from, keeping transformation F
[T,F] = reduction(S);
%save('Tshift.mat','T');
[m,n] = size(T);
H = eye(n,n);
D = eye(n,n);

%Keep track of iterations
q = 0;

%Convergence criterion
while (norm(S*(F*H)-(F*H)*D) > 10^(-8))
    
    %Iterate until bottom subdiagonal entry is nearly zero, 
    %meaning bottom diagonal entry is nearly an eigenvalue
    while (n > 1 && abs(T(n,n-1)) > eps*(abs(T(n,n))+abs(T(n-1,n-1)))) 
        I = eye(n);
        l = T(n,n);
        T = T - l*I;
        [Q,A] = qrIterateTridiag(T);
        T = A + l*I;
        disp(T);
        M = eye(m);
        M(1:n,1:n) = Q;
        H = H*M;
        q = q+1;
    end
    
    %Deflate once we find an approximate eigenvalue
    if(n > 2)
        D(n,n) = T(n,n);
        n = n-1;
        T = T(1:n,1:n);
        disp(norm(S*F*H-F*H*D));
        disp(n);
        disp(q);
    end
    
    %If T is now 2x2, find eigenvalue by hand
    if n == 2
        [V,B] = eigs(T);
        D(1,1) = B(1,1);
        D(2,2) = B(2,2);
        M = eye(m,m);
        M(1:2,1) = V(1:2,2);
        M(1:2,2) = V(1:2,1);
        disp(D);
        H = H*M;
        n = 1;
    end
    
    %Form matrix of eigenvectors
    P = F*H;
    
    %Save results
    %save('Pshift.mat','P');
    %save('Dshift.mat','D');
    %save('Hshift.mat','H');
    disp(['Norm SP-PD: ' num2str(norm(S*P-P*D))]);
    disp(['Number of QR with shift iterations: ' num2str(q)]);
end
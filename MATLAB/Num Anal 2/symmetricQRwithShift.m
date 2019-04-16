%Joshua Enxing
%University of Connecticut
%MA5511

%Does implicit-Q algorithm for a real symmetric matrix S,
%first reducing it to a symmetric tridiagonal matrix T
function symmetricQRwithShift(S)
tic
[T,F] = reduction(S);
save('Tshift.mat','T');
[m,n] = size(T);
H = eye(n);
D = eye(n);
q = 0;
while (n > 1)
    if n > 2
        while ((abs(T(n,n-1)) > eps*(abs(T(n,n))+ abs(T(n-1,n-1)))))
            %Get first column of T to do implcit Q
            A = T(1:2,1);
            l = T(n,n);
            A(1) = A(1) - l;
            
            %Find values for givens matrix based on 1st column of 
            %shifted matrix
            if (((A(1))^2 + (A(2))^2) ~= 0)
                c = A(1)/sqrt(A(1)^2 + A(2)^2);
                s = A(2)/sqrt(A(1)^2 + A(2)^2);
            else
                c = 1;
                s = 0;
            end
            
            %Make givens matrix 
            G(1,1) = c;
            G(1,2) = s;
            G(2,1) = -s;
            G(2,2) = c;
            
            %Left multiplication
            T(1:2,1:n) = G*T(1:2,1:n);
            
            G(2,1) = s;
            G(1,2) = -s;
            
            %Right mult by G'
            T(1:n,1:2) = T(1:n,1:2)*G;
            
            %Update H
            H(1:m,1:2) = H(1:m,1:2)*G;
            
            %Chase the bulge
            for j=2:n-1
                
                %Find values for givens matrix to get rid of bulge
                if (((T(j,j-1))^2 + (T(j+1,j-1))^2) ~= 0)
                    c = T(j,j-1)/sqrt((T(j,j-1))^2 + (T(j+1,j-1))^2);
                    s = T(j+1,j-1)/sqrt((T(j,j-1))^2 + (T(j+1,j-1))^2);
                else
                    c = 1;
                    s = 0;
                end
          
                %Form givens matrix
                G(1,1) = c;
                G(1,2) = s;
                G(2,1) = -s;
                G(2,2) = c;
                
                %left mult
                T(j:j+1,1:n) = G*T(j:j+1,1:n);
                
                G(2,1) = s;
                G(1,2) = -s;
                %right mult
                T(1:n,j:j+1) = T(1:n,j:j+1)*G;
                
                %T = T*G;
                H(1:m,j:j+1) = H(1:m,j:j+1)*G;
            end
            
            q = q+1;
        end
        
    end
    
    %Deflate nxn matrix to (n-1)x(n-1)
    if(n > 2)
        D(n,n) = T(n,n);
        n = n-1;
        T = T(1:n,1:n);
    end
    
    %Diagonalize the final 2x2
    if n == 2
        while (abs(T(1,2)) > (10^(-8))*(abs(T(2,2))+abs(T(2,1))))
            [Q,T] = qrIterateTridiag(T);
            H(1:m,1:2) = H(1:m,1:2)*Q;
        end
        D(1,1) = T(1,1);
        D(2,2) = T(2,2);
        n = 1;
    end
end
save('Pshift.mat','P');
save('Dshift.mat','D');
save('Hshift.mat','H');
P = F*H;
disp(['Norm SP-PD: ' num2str(norm(S*P-P*D))]);
disp(['Number of implicit-Q iterations: ' num2str(q)]);
disp(['Time to run: ' num2str(toc)]);
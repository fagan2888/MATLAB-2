%Joshua Enxing
%University of Connecticut
%MA5511

%Does standard QR algorithm on a real symmetric matrix S, 
%first reducing it to a symmetric tridiagonal matrix T
function symmetricQR(S)
tic
%Reduce to symmetric tridiagonal
[T,F] = reduction(S);
save('Tnoshift.mat','T');
[m,n] = size(T);

H = eye(n);
D = eye(n);

%Form matrix of eigenvectors. H will transform T to a nearly-diagonal
%matrix, since each iteration will result in a symmetric matrix. F*H will
%transform eigenvectors of this diagonal matrix (standard basis vectors)
%into eigenvectors of S
P = F*H;
q = 0;
while (norm(S*P - P*D) > 10^(-8))    
    %Do a QR iteration. Factors T into T = QR then sets
    %T = Q'*T*Q, without explicitly forming R
    [Q,T] = qrIterateTridiag(T);
    q = q+1;
    %H here is the composition of all the QR iterations so we 
    % can back-transform to get eigenvectors
    H = H*Q;
    for i=1:n
        D(i,i) = T(i,i);
    end
    P = F*H;
end

%Save results
save('Pnoshift.mat','P');
save('Dnoshift.mat','D');
save('Hnoshift.mat','H');
disp(['Norm SP-PD: ' num2str(norm(S*P-P*D))]);
disp(['Number of QR iterations: ' num2str(q)]);
disp(['Time to run: ' num2str(toc)]);
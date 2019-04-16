%Weighted Jacobi method for a tridiagonal matrix A
function x = weighted_jacobi_tridiag(A,x_guess,b_vect,w_weight)

x = x_guess;
b = b_vect;
w = w_weight;
n = length(x);

temp = x;

for i=1:n
    temp(i) = b(i)/A(i,i);
    %Due to the tridiagonal structure of A
    for j=i-1:i+1
        if j~=0 && j~=n+1 && j~=i
            temp(i) = temp(i) - (A(i,j)*x(j))/A(i,i);
        end
    end
    temp(i) = (1-w)*x(i) + w*temp(i);
end

x = temp;
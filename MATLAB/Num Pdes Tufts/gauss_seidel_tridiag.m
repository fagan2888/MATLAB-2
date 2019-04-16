%Gauss_Seidel method for a tridiagonal matrix A
function x = gauss_seidel_tridiag(A,x_guess,b_vect)

x = x_guess;
b = b_vect;
n = length(x);

for i=1:n
    x(i) = b(i)/A(i,i);
    %Due to the tridiagonal structure of A
    for j=i-1:i+1
        if j~=0 && j~=n+1 && j~=i
            x(i) = x(i) - (A(i,j)*x(j))/A(i,i);
        end
    end
end


function [A,b] = make_matrix_A_and_vector_b(n)

A = zeros(n);

for i=1:n
    for j=1:n
        A(i,j) = 1/(i+j-1);
    end
end

b = ones(n,1);

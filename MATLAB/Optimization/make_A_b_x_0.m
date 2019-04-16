%Joshua Enxing
%Tufts University
%MA150

%Produces matrix A, vector b, and vector x_0 for 
%Problem 1 in HW 2 based on input dimension n
function [A,b,x] = make_A_b_x_0(n)

A = zeros(n);

for i=1:n
    for j=1:n
        A(i,j) = 1/(i+j-1);
    end
end

b = ones(n,1);
x = zeros(n,1);

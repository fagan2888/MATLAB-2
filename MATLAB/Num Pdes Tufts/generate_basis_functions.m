function thetas = generate_basis_functions(partition)
y = partition;
n = length(y);

thetas = sym(zeros(n,1));

thetas(1) = @(x)heaviside(y(2)-x).*((x-y(2))/(y(1)-y(2)));

for i=1:n-2
    thetas(i+1) = @(x)heaviside(x-y(i)).*((x-y(i))/(y(i+1)-y(i))).*heaviside(y(i+1)-x)+heaviside(x-y(i+1)).*((x-y(i+2))/(y(i+1)-y(i+2))).*heaviside(y(i+2)-x);
end

thetas(n) = @(x)heaviside(x-y(n-1)).*((x-y(n-1))/(y(n)-y(n-1)));

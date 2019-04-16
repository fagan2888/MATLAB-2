function x = gauss_seidel(A,x,b)

n = length(x);
m = sqrt(n);
% for i = 1:n
%     x(i) = (1/A(i,i))*(b(i) - A(i,:)*x + A(i,i)*x(i));
% end

x(1) = (1/A(1,1))*(b(1) - A(1,2)*x(2) - A(1,m+1)*x(m+1));
for i = 2:n-1
    if i > m && i < n-m+1
        x(i) = (1/A(i,i))*(b(i) - A(i,i-m)*x(i-m) -  A(i,i-1)*x(i-1) - A(i,i+1)*x(i+1) - A(i,i+m)*x(i+m));
    elseif i <= m
        x(i) = (1/A(i,i))*(b(i) - A(i,i-1)*x(i-1) - A(i,i+1)*x(i+1) - A(i,i+m)*x(i+m));
    elseif i >= n-m+1
        x(i) = (1/A(i,i))*(b(i) - A(i,i-m)*x(i-m) -  A(i,i-1)*x(i-1) - A(i,i+1)*x(i+1));
    end
end
x(n) = (1/A(n,n))*(b(n) - A(n,n-m)*x(n-m) -  A(n,n-1)*x(n-1));

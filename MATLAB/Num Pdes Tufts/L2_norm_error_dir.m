function e = L2_norm_error_tentc2(u,u_true,partition)

e = 0;
y = partition;
n = length(y);

%on [y(1),y(2)], our approx to u is u(2)*phi_2 since u(0) = 0
f = @(x)(u(2)*((x-y(1))/(y(2)-y(1)))-u_true(x)).^2;
e = e + integral(f,y(1),y(2));

%Loop over all eles except first and last
for i=2:n-2
    %on [y(i),y(i+1)], our approx to u is u(i)*phi_i + u(i+1)*phi_(i+1)
    f = @(x)(u(i)*((x-y(i+1))/(y(i)-y(i+1)))+u(i+1)*((x-y(i))/(y(i+1)-y(i)))-u_true(x)).^2;
    e = e + integral(f,y(i),y(i+1));
end

%on [y(n-1),y(n)], our approx to u is u(n-1)*phi_(n-1) since u(1) = 0
f = @(x)(u(n-1)*((x-y(n))/(y(n-1)-y(n)))-u_true(x)).^2;
e = e + integral(f,y(n-1),y(n));


e = sqrt(e);
function e = L2_norm_error(u,u_true,partition)

e = 0;
y = partition;
n = length(y);

%Loop over all eles
for i=1:n-1
    %on [y(i),y(i+1)], our approx to u is u(i)*phi_i + u(i+1)*phi_(i+1)
    f = @(x)(u(i)*((x-y(i+1))/(y(i)-y(i+1)))+u(i+1)*((x-y(i))/(y(i+1)-y(i)))-u_true(x)).^2;
    e = e + integral(f,y(i),y(i+1));
end

e = sqrt(e);
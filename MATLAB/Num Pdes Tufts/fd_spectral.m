function error = fd_spectral(N)

x = linspace(-1,1,N);
dx = x(2) - x(1);
u = zeros(N,1);

num_steps = 2*ceil(1/dx);
dt = 1/num_steps;
lambda = dt/dx;


Q = full(gallery('tridiag',N,0,1-lambda,lambda));

disp(dx);
disp(dt);

for i=1:N
    u(i) = sin(pi*x(i));
end

for i=1:num_steps
    u = Q*u;
    u(N) = sin(pi*(x(N)+i*dt));
end

sum = 0;
for i=1:N
    sum = sum + (u(i)-sin(pi*(x(i)+1)))^2;
end
error = sqrt(sum);
plot(x,u)

    

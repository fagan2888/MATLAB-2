function [y,u,error] = chebyshev_spectral_dirichlet(init,N,t_f,dt)

syms x;
y = make_collocation_points(N);
u_true = zeros(N,1);
D = poldif(y,2);
D2 = D(:,:,1);

f = @(t) sin(pi*(1+t));
u_fun = @(x,t) sin(pi*(x+t));
u = zeros(N,1);

for i=1:N
    u(i) = init(y(i));
    u_true(i) = u_fun(y(i),t_f);
end

while norm(D2) >= 2/dt
    dt = dt/2;
end
disp(dt);

num_steps = t_f/dt;
for i=1:num_steps
    u = (eye(N)-(dt/2)*D2)\((eye(N)+(dt/2)*D2)*u);
    u(1) = f(dt*i);
end

sum = 0;
for i=1:N
    sum = sum + abs(u(i)-u_true(i))^2;
end
error = sqrt(sum);

plot(y,u);

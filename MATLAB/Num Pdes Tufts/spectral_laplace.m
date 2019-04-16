function u = spectral_laplace(N,f)

F = zeros(N,1);
A = zeros(N,1);
u = zeros(N,1);

u(1) = integral(f,0,2*pi)/(2*pi);
for i=2:N
    fun = @(x) sin((i-1)*x);
    g = @(x) f(x).*fun(x);
    F(i) = integral(g,0,2*pi);
    
    A(i) = ((i-1)^2+1)*pi;
    
    u(i) = F(i)/A(i);
end


u_fun = @(x) u(1);
for i=2:N
    u_fun = @(x) u_fun(x) + u(i)*sin(i*x);
end

u_fun = @(x) u_fun(x) + u(N);

fplot(u_fun,[0 2*pi]);






    
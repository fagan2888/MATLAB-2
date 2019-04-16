function u_fun = FE_neumann(N,f)

F = zeros(N,1);
A = zeros(N);
u = zeros(N,1);
y = linspace(-1,1,N);
hat = generate_basis_functions(y);
h = y(2) - y(1);

syms x;

for i=2:N-1
    F(i) = int(f*hat(i),x,y(i)-h,y(i)+h);
end
F(1) = int(f*hat(1),x,-1,h-1);
F(N) = int(f*hat(N),x,1-h,1);

A(1,1) = int(hat(1)*hat(1) + (1/h^2),x,y(1),y(2));
A(1,2) = int(hat(1)*hat(2) - (1/h^2),x,y(1),y(2));
for i=2:N-1
    A(i,i-1) = A(i,i-1) + int(hat(i)*hat(i-1) - (1/h^2),x,y(i-1),y(i));
    A(i,i) = A(i,i) + int(hat(i)*hat(i) + (1/h^2),x,y(i-1),y(i+1));
    A(i,i+1) = A(i,i+1) + int(hat(i)*hat(i+1) - (1/h^2),x,y(i),y(i+1));
end
A(N,N-1) = A(N,N-1) + int(hat(N)*hat(N-1) - (1/h^2),x,y(N-1),y(N));
A(N,N) = A(N,N) + int(hat(N)*hat(N) + (1/h^2),x,y(N-1),y(N));

u = A\F;

u_fun = u(1)*hat(1);
for i=2:N
    u_fun = u_fun + u(i)*hat(i);
end
 
fplot(u_fun,[-1 1])
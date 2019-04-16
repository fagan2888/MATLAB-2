function u_fun = spectral_legendre(N,f)

F = zeros(N,1);
A = zeros(N);
u = zeros(N,1);
L = sym(zeros(N,1));
D = sym(zeros(N,1));

syms x;

for i=1:N
    L(i) = legendreP(i-1,x);
    D(i) = diff(L(i),x);
    F(i) = int(f*L(i),x,-1,1);
end

for i=1:N
    for j=1:N
        A(i,j) = int(D(i)*D(j)+L(i)*L(j),x,-1,1);
    end
end

u = A\F;

u_fun = u(1)*L(1);
for i=2:N
    u_fun = u_fun + u(i)*L(i);
end

fplot(u_fun,[-1 1])




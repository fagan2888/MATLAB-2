function theta = generate_thetas(z, n)
syms x;
theta = zeros(n+2,1);

theta(1) = (x-z(2))/(0-z(2));
theta(n+2) = (x-z(n+1))/(1-z(n+1));

for i=2:n+1
    theta(i) = heaviside(x(i)-x)*((x-x(i-1))/(x(i)-x(i-1))) + heaviside(x-x(i))*((x-x(i+1))/(x(i)-x(i+1)));
end

ezplot(theta(1));
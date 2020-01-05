function g = grad_f(x)

n = length(x);
g = zeros(n,1);

for i=1:n
    if mod(i,2)==0
        g(i) = 200*(x(i)-x(i-1)^2);
    else
        g(i) = -400*x(i)*(x(i+1)-x(i)^2)-2*(1-x(i));
    end
end
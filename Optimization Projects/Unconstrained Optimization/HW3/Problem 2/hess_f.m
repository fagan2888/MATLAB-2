function h = hess_f(x)

n = length(x);
h1 = zeros(n,1);
h2 = zeros(n-1,1);

for i=1:n-1
    if mod(i,2)==0
        h1(i) = 200;
        h2(i) = 0;
    else
        h1(i) = -400*(x(i+1)-3*x(i)^2)+2;
        h2(i) = -400*x(i);
    end
end
h1(n) = 200;

h = full(gallery('tridiag',h2,h1,h2));
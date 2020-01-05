function f = extended_rosenbrock(x)

n = length(x);
f = 0;

for i=1:0.5*n
    f = f + 100*(x(2*i)-x(2*i-1)^2)^2+(1-x(2*i-1))^2;
end